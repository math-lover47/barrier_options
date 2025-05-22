classdef FEMSolver < handle
    % FEMSolver - Finite Element Solver for Black-Scholes PDE (European Option)
    % Supports linear (P1) and quadratic (P2) finite elements
    % Handles time-dependent boundary conditions correctly during time stepping

    properties
        ne              % Number of elements
        nodes           % Asset price grid (spatial nodes)
        dof             % Degrees of freedom
        M               % Mass matrix
        K               % Diffusion matrix (with S^2 scaling)
        C               % Convection matrix (with S scaling)
        R               % Reaction matrix (r * mass)
        theta           % Time-stepping parameter
        u0              % Solution vector (option values at nodes)
        x_end           % Maximum asset price
        connectivity    % Element connectivity
        element_type    % 'P1' or 'P2'

        % Black-Scholes parameters
        K_strike
        sigma
        r
        optionType
        T
        current_time    % Track current time in backward stepping
    end

    methods
        function obj = FEMSolver(ne, element_type, x_end, theta, K_strike, sigma, r, T, optionType)
            obj.ne = ne;
            obj.element_type = element_type;
            obj.x_end = x_end;
            obj.theta = theta;
            obj.K_strike = K_strike;
            obj.sigma = sigma;
            obj.r = r;
            obj.T = T;
            obj.optionType = lower(optionType);
            obj.current_time = T;

            obj.initialize_nodes();
            obj.assemble_system();
            obj.set_payoff_initial_condition();
        end

        function initialize_nodes(obj)
            if strcmp(obj.element_type,'P1')
                obj.nodes = linspace(0,obj.x_end,obj.ne+1)';
                obj.dof = obj.ne+1;
                obj.connectivity = [(1:obj.ne)' (2:obj.ne+1)'];
            elseif strcmp(obj.element_type,'P2')
                obj.nodes = linspace(0,obj.x_end,2*obj.ne+1)';
                obj.dof = 2*obj.ne+1;
                obj.connectivity = zeros(obj.ne,3);
                for e=1:obj.ne
                    obj.connectivity(e,:) = [2*e-1,2*e+1,2*e]; % left, right, middle
                end
            else
                error('Unsupported element type.');
            end
        end

        function assemble_system(obj)
            obj.M = sparse(obj.dof,obj.dof);
            obj.K = sparse(obj.dof,obj.dof);
            obj.C = sparse(obj.dof,obj.dof);
            obj.R = sparse(obj.dof,obj.dof);
            if strcmp(obj.element_type,'P1')
                obj.assemble_p1_system();
            else
                obj.assemble_p2_system();
            end
        end

        function assemble_p1_system(obj)
            [pts,wts] = obj.get_quadrature(3);
            for e=1:obj.ne
                idx = obj.connectivity(e,:);
                x_e = obj.nodes(idx);
                h = x_e(2)-x_e(1);
                M_e = zeros(2); K_e = zeros(2); C_e = zeros(2); R_e = zeros(2);
                for q=1:length(wts)
                    xi = pts(q); w = wts(q);
                    x_phys = x_e(1) + (xi+1)*h/2;
                    J = h/2;
                    N = [(1-xi)/2, (1+xi)/2];
                    dNdx = [-1/2, 1/2]/J;
                    if x_phys < 1e-12
                        D = 0; Vc = 0;
                    else
                        D = 0.5*obj.sigma^2*x_phys^2;
                        Vc = obj.r*x_phys;
                    end
                    % Mass
                    M_e = M_e + w*J*(N' * N);
                    % Reaction (r * mass)
                    R_e = R_e + w*J*obj.r*(N' * N);
                    % Diffusion
                    K_e = K_e + w*J*D*(dNdx' * dNdx);
                    % Convection
                    C_e = C_e + w*J*Vc*(N' * dNdx);
                end
                obj.M(idx,idx) = obj.M(idx,idx) + M_e;
                obj.R(idx,idx) = obj.R(idx,idx) + R_e;
                obj.K(idx,idx) = obj.K(idx,idx) + K_e;
                obj.C(idx,idx) = obj.C(idx,idx) + C_e;
            end
        end

        function assemble_p2_system(obj)
            [pts,wts] = obj.get_quadrature(4);
            for e=1:obj.ne
                idx = obj.connectivity(e,:);
                x_e = obj.nodes(idx);
                h = x_e(2) - x_e(1);
                M_e=zeros(3); K_e=zeros(3); C_e=zeros(3); R_e=zeros(3);
                for q=1:length(wts)
                    xi = pts(q); w = wts(q);
                    x_phys = x_e(1) + (xi+1)*h/2;
                    J = h/2;
                    N = [xi*(xi-1)/2, xi*(xi+1)/2, 1-xi^2];
                    dNdx = [(2*xi-1)/2, (2*xi+1)/2, -2*xi] / J;
                    if x_phys < 1e-12
                        D = 0; Vc = 0;
                    else
                        D = 0.5*obj.sigma^2*x_phys^2;
                        Vc = obj.r*x_phys;
                    end
                    M_e = M_e + w*J*(N' * N);
                    R_e = R_e + w*J*obj.r*(N' * N);
                    K_e = K_e + w*J*D*(dNdx' * dNdx);
                    C_e = C_e + w*J*Vc*(N' * dNdx);
                end
                obj.M(idx,idx) = obj.M(idx,idx) + M_e;
                obj.R(idx,idx) = obj.R(idx,idx) + R_e;
                obj.K(idx,idx) = obj.K(idx,idx) + K_e;
                obj.C(idx,idx) = obj.C(idx,idx) + C_e;
            end
        end

        function set_payoff_initial_condition(obj)
            obj.u0 = zeros(obj.dof,1);
            switch obj.optionType
                case 'call'
                    obj.u0 = max(obj.nodes - obj.K_strike,0);
                case 'put'
                    obj.u0 = max(obj.K_strike - obj.nodes,0);
            end
        end

        function [A,b] = apply_boundary_conditions(obj,A,b,time)
            tau = obj.T - max(time,0);
            switch obj.optionType
                case 'call'
                    bcL = 0;
                    bcR = max(0, obj.nodes(end) - obj.K_strike*exp(-obj.r*tau));
                case 'put'
                    bcL = obj.K_strike*exp(-obj.r*tau);
                    bcR = 0;
            end
            A(1,:) = 0; A(1,1) = 1; b(1) = bcL;
            A(end,:) = 0; A(end,end) = 1; b(end) = bcR;
        end

        function u_next = step(obj,dt)
            obj.current_time = obj.current_time - dt;
            L = obj.K + obj.C - obj.R;  % Combine with correct reaction sign
            A = obj.M + obj.theta*dt*L;
            b = (obj.M - (1-obj.theta)*dt*L) * obj.u0;
            [A,b] = obj.apply_boundary_conditions(A,b,obj.current_time);
            obj.u0 = A \ b;
            u_next = obj.u0;
        end

        function [pts,wts] = get_quadrature(~,n)
            switch n
                case 1, pts = 0; wts = 2;
                case 2, pts = [-1/sqrt(3), 1/sqrt(3)]; wts = [1,1];
                case 3, pts = [-sqrt(3/5), 0, sqrt(3/5)]; wts = [5/9,8/9,5/9];
                case 4
                    pts = [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7), ...
                           sqrt((3-2*sqrt(6/5))/7),  sqrt((3+2*sqrt(6/5))/7)];
                    wts = [(18-sqrt(30))/36, (18+sqrt(30))/36, ...
                           (18+sqrt(30))/36, (18-sqrt(30))/36];
                otherwise
                    error('Quadrature rule for n=%d not implemented', n);
            end
        end
    end
end

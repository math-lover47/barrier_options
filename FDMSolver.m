classdef FDMSolver < handle
    % FDMSolver - Finite Difference Solver for Black-Scholes PDE (European Option)
    % Solves: dV/dt + 0.5*sigma^2*S^2*d2V/dS2 + r*S*dV/dS - r*V = 0
    % Using Crank-Nicolson time stepping
    
    properties
        nx          % Number of spatial points
        dx          % Spatial step size
        xvec        % Spatial grid (asset prices)
        theta       % Theta parameter for Crank-Nicolson (0.5)
        u0          % Current solution vector
        x_end       % Spatial domain end (max asset price)
        K           % Strike price (for boundary conditions)
        r           % Risk-free rate (for boundary conditions)
        optionType  % 'call' or 'put' (for boundary conditions)
    end
    
    methods
        function obj = FDMSolver(ne, x_end, theta)
            obj.nx = ne + 1;
            obj.dx = x_end / ne;
            obj.xvec = linspace(0, x_end, obj.nx)';
            obj.theta = theta;
            obj.x_end = x_end;
        end
        
        function setInitialCondition(obj, ic_fun, K, r, optionType)
            % Set initial condition at maturity (payoff function)
            % Also store parameters needed for boundary conditions
            obj.u0 = ic_fun(obj.xvec);
            obj.K = K;
            obj.r = r;
            obj.optionType = lower(optionType);
        end
        
        function u_next = step(obj, t_current, dt, sigma, r)
            % One Crank-Nicolson step backward in time for Black-Scholes PDE
            % The PDE is: dV/dt + 0.5*sigma^2*S^2*d2V/dS2 + r*S*dV/dS - r*V = 0
            % We solve: dV/dt = -0.5*sigma^2*S^2*d2V/dS2 - r*S*dV/dS + r*V
            
            nx_interior = obj.nx - 2;
            S = obj.xvec(2:end-1); % Interior grid points (asset prices)
            u_interior = obj.u0(2:end-1); % Interior solution values
            
            % Variable coefficients at interior points
            D = 0.5 * sigma^2 * S.^2;  % Diffusion coefficients
            C = r * S;                  % Convection coefficients
            
            % Build finite difference matrices for interior points
            % Second derivative matrix (d2V/dS2)
            K_matrix = obj.buildSecondDerivativeMatrix(nx_interior, obj.dx);
            
            % First derivative matrix (dV/dS)
            C_matrix = obj.buildFirstDerivativeMatrix(nx_interior, obj.dx);
            
            % Apply variable coefficients
            D_diag = spdiags(D, 0, nx_interior, nx_interior);
            C_diag = spdiags(C, 0, nx_interior, nx_interior);
            
            % Assemble the spatial operator: L = D*K_matrix + C_diag*C_matrix - r*I
            L = D_diag * K_matrix + C_diag * C_matrix - r * speye(nx_interior);
            
            % Identity matrix
            I = speye(nx_interior);
            
            % Crank-Nicolson matrices
            % (I - theta*dt*L) * u^{n+1} = (I + (1-theta)*dt*L) * u^n + boundary terms
            A = I - obj.theta * dt * L;
            B = I + (1 - obj.theta) * dt * L;
            
            % Right hand side
            rhs = B * u_interior;
            
            % Add boundary condition contributions
            tau = t_current; % Time to maturity
            [bc_left, bc_right] = obj.getBoundaryConditions(tau);
            
            % Boundary contributions to the right-hand side
            % Left boundary (S = 0)
            if nx_interior > 0
                % Contribution from left boundary to first interior point
                D_0 = 0.5 * sigma^2 * 0^2; % = 0
                C_0 = r * 0; % = 0
                
                % Second derivative contribution: D_0 * (bc_left - 2*u_1 + u_2) / dx^2
                % First derivative contribution: C_0 * (u_2 - bc_left) / (2*dx)
                bc_contrib_left = 0; % Both D_0 and C_0 are zero at S=0
                
                rhs(1) = rhs(1) + bc_contrib_left;
            end
            
            % Right boundary (S = S_max)
            if nx_interior > 0
                S_max = obj.xvec(end);
                D_max = 0.5 * sigma^2 * S_max^2;
                C_max = r * S_max;
                
                % Second derivative: D_max * (u_{n-1} - 2*u_n + bc_right) / dx^2
                % First derivative: C_max * (bc_right - u_{n-1}) / (2*dx)
                
                % Theta scheme contributions
                bc_contrib_right_old = (1 - obj.theta) * dt * ...
                    (D_max * bc_right / (obj.dx^2) + C_max * bc_right / (2 * obj.dx));
                bc_contrib_right_new = obj.theta * dt * ...
                    (D_max * bc_right / (obj.dx^2) + C_max * bc_right / (2 * obj.dx));
                
                rhs(end) = rhs(end) + bc_contrib_right_old;
                
                % Modify matrix A for new time step boundary contribution
                if strcmpi(obj.optionType, 'call')
                    A(end, end) = A(end, end) - bc_contrib_right_new / bc_right;
                end
            end
            
            % Solve the linear system
            u_next_interior = A \ rhs;
            
            % Assemble full solution with boundary conditions
            u_next = zeros(obj.nx, 1);
            u_next(1) = bc_left;
            u_next(2:end-1) = u_next_interior;
            u_next(end) = bc_right;
            
            % Update stored solution
            obj.u0 = u_next;
        end
        
        function [bc_left, bc_right] = getBoundaryConditions(obj, tau)
            % Get boundary conditions at time tau (time to maturity)
            % For European options:
            % Call: V(0,t) = 0, V(S_max,t) ≈ S_max - K*exp(-r*tau)
            % Put: V(0,t) = K*exp(-r*tau), V(S_max,t) = 0
            
            if strcmp(obj.optionType, 'call')
                bc_left = 0; % Call option worth 0 when S = 0
                bc_right = max(0, obj.x_end - obj.K * exp(-obj.r * tau)); % Deep ITM call
            else % put option
                bc_left = obj.K * exp(-obj.r * tau); % Put option when S = 0
                bc_right = 0; % Put option worth 0 when S >> K
            end
        end
        
        function A = buildSecondDerivativeMatrix(~, N, dx)
            % Build second derivative matrix using central differences
            % d2u/dx2 ≈ (u_{i-1} - 2*u_i + u_{i+1}) / dx^2
            e = ones(N, 1);
            A = spdiags([e -2*e e], -1:1, N, N) / (dx^2);
        end
        
        function A = buildFirstDerivativeMatrix(~, N, dx)
            % Build first derivative matrix using central differences
            % du/dx ≈ (u_{i+1} - u_{i-1}) / (2*dx)
            e = ones(N, 1);
            A = spdiags([-e zeros(N,1) e], -1:1, N, N) / (2 * dx);
        end
    end
end
classdef option_new < handle
    % Class for pricing options under the Black-Scholes model.
    %
    % Assumptions:
    %   - Constant risk-free rate
    %   - Constant volatility
    %   - Constant continuous dividend yield

    properties
        r              % Risk-free rate
        q              % Continuous dividend yield
        spot_price     % Current price of the underlying asset
        strike         % Strike price
        sig          % Volatility
        maturity       % Time to maturity (in years)
        option_type    % Type of option: 'VANILLA', 'BARRIER', etc.
        exercise_type  % 'EUROPEAN' or 'AMERICAN'
        position       % 'CALL' or 'PUT'
        barrier        % Single barrier level
        lower_barrier  % Lower barrier (double barrier options)
        upper_barrier  % Upper barrier (double barrier options)
        util           % Utility functions object
        theta          % FDM theta-scheme parameter
    end

    methods
        function obj = option_new(r, q, spot_price, strike, sig, t, ...
                                  option_type, exercise_type, position, ...
                                  t_switch, theta, varargin)
            % Constructor for option_new class

            obj.util = utils();
            obj.r = r;
            obj.q = q;
            obj.spot_price = spot_price;
            obj.strike = strike;
            obj.sig = sig;
            obj.theta = theta;

            % Handle time interpretation
            if nargin < 10 || isempty(t_switch)
                t_switch = 0;
            end

            switch t_switch
                case 0
                    obj.maturity = t / 252;   % Trading days
                case 1
                    obj.maturity = t / 365;   % Calendar days
                case 2
                    obj.maturity = t;         % In years
                otherwise
                    error('Invalid t_switch value');
            end

            % Validate and set option type
            valid_types = {'VANILLA', 'ASIAN', 'LOOKBACK', ...
                'UP-AND-IN-BARRIER', 'UP-AND-OUT-BARRIER', ...
                'DOWN-AND-IN-BARRIER', 'DOWN-AND-OUT-BARRIER', ...
                'KNOCK-OUT-DOUBLE-BARRIER', 'KNOCK-IN-DOUBLE-BARRIER'};
            option_type = upper(option_type);

            if ismember(option_type, valid_types)
                obj.option_type = option_type;
            else
                error('Unsupported option type: %s', option_type);
            end

            % Set exercise type
            exercise_type = upper(exercise_type);
            if ismember(exercise_type, {'EUROPEAN', 'AMERICAN'})
                obj.exercise_type = exercise_type;
            else
                error('Invalid exercise type: %s', exercise_type);
            end

            % Set position
            position = upper(position);
            if ismember(position, {'CALL', 'PUT'})
                obj.position = position;
            else
                error('Invalid position type: %s', position);
            end

            % Parse optional parameters
            p = inputParser;
            addParameter(p, 'barrier', []);
            addParameter(p, 'lower_barrier', []);
            addParameter(p, 'upper_barrier', []);
            parse(p, varargin{:});

            % Validate barriers
            switch obj.option_type
                case {'UP-AND-IN-BARRIER', 'UP-AND-OUT-BARRIER'}
                    if isempty(p.Results.barrier)
                        error('Missing barrier value');
                    elseif p.Results.barrier <= obj.spot_price
                        error('Barrier must be greater than spot price');
                    else
                        obj.barrier = p.Results.barrier;
                    end

                case {'DOWN-AND-IN-BARRIER', 'DOWN-AND-OUT-BARRIER'}
                    if isempty(p.Results.barrier)
                        error('Missing barrier value');
                    elseif p.Results.barrier >= obj.spot_price
                        error('Barrier must be less than spot price');
                    else
                        obj.barrier = p.Results.barrier;
                    end

                case {'KNOCK-IN-DOUBLE-BARRIER', 'KNOCK-OUT-DOUBLE-BARRIER'}
                    lb = p.Results.lower_barrier;
                    ub = p.Results.upper_barrier;
                    if isempty(lb) || isempty(ub)
                        error('Missing lower/upper barrier values');
                    elseif ub <= lb
                        error('Upper barrier must be greater than lower barrier');
                    elseif obj.spot_price >= ub || obj.spot_price <= lb
                        error('Spot price must lie between the barriers');
                    else
                        obj.lower_barrier = lb;
                        obj.upper_barrier = ub;
                    end
            end
        end

        function price = bs_call(obj)
            % Price a European call option using the Black-Scholes formula
            d1 = (log(obj.spot_price / obj.strike) + ...
                (obj.r - obj.q + 0.5 * obj.sig^2) * obj.maturity) / ...
                (obj.sig * sqrt(obj.maturity));
            d2 = d1 - obj.sig * sqrt(obj.maturity);

            forward = obj.spot_price * exp((obj.r - obj.q) * obj.maturity);
            price = exp(-obj.r * obj.maturity) * ...
                (forward * option.norm_cdf(d1) - obj.strike * option.norm_cdf(d2));
        end

        function price = bs_put(obj)
            % Price a European put option using the Black-Scholes formula
            d1 = (log(obj.spot_price / obj.strike) + ...
                (obj.r - obj.q + 0.5 * obj.sig^2) * obj.maturity) / ...
                (obj.sig * sqrt(obj.maturity));
            d2 = d1 - obj.sig * sqrt(obj.maturity);

            forward = obj.spot_price * exp((obj.r - obj.q) * obj.maturity);
            price = exp(-obj.r * obj.maturity) * ...
                (obj.strike * option.norm_cdf(-d2) - forward * option.norm_cdf(-d1));
        end

        function v0 = fdm_vanilla(obj, Ns, Nt, m, theta)
            % FDM pricing for vanilla option using theta scheme

            Nt = Nt * m;  % Extend time grid
            mu = obj.r - obj.q;
            range = 3 * obj.sig * sqrt(obj.maturity);

            Smax = obj.spot_price * exp((mu - 0.5 * obj.sig^2) * obj.maturity + range);
            Smin = 0;

            dt = obj.maturity / Nt;
            ds = (Smax - Smin) / Ns;
            sGrid = linspace(Smin, Smax, Ns + 1);

            % Terminal payoff
            if strcmpi(obj.position, 'CALL')
                V = max(sGrid - obj.strike, 0);
                f0 = zeros(Nt + 1, 1);
                fNs = Smax * exp(-obj.q * linspace(0, obj.maturity, Nt + 1)') - ...
                      obj.strike * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)');
            elseif strcmpi(obj.position, 'PUT')
                V = max(obj.strike - sGrid, 0);
                f0 = obj.strike * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)');
                fNs = zeros(Nt + 1, 1);
            else
                error('Invalid position type');
            end

            i = (1:Ns-1)';
            alpha = 0.5 * dt * (obj.sig^2 * i.^2 - mu * i);
            beta  = 1 + dt * (obj.sig^2 * i.^2 + obj.r);
            gamma = 0.5 * dt * (obj.sig^2 * i.^2 + mu * i);

            a_A = -theta * alpha;
            b_A = beta - 2 * theta * dt * obj.sig^2 * i.^2;
            c_A = -theta * gamma;

            a_B = (1 - theta) * alpha;
            b_B = beta - 2 * (1 - theta) * dt * obj.sig^2 * i.^2;
            c_B = (1 - theta) * gamma;

            V_curr = V(2:Ns);  % Interior values only

            for k = Nt:-1:1
                rhs = zeros(Ns - 1, 1);
                for j = 2:Ns-2
                    rhs(j) = a_B(j) * V_curr(j - 1) + ...
                             b_B(j) * V_curr(j) + ...
                             c_B(j) * V_curr(j + 1);
                end
                % Boundaries
                rhs(1)     = rhs(1)     + a_B(1)   * V_curr(1)   + c_B(1)   * V_curr(2)   - a_A(1)   * f0(k);
                rhs(end)   = rhs(end)   + a_B(end) * V_curr(end-1) + c_B(end) * V_curr(end) - c_A(end) * fNs(k);

                ab = obj.util.compress_tridiag(a_A, b_A, c_A);
                V_next = obj.util.thomas_algorithm(ab, rhs);

                V_curr = V_next;
            end

            % Interpolation for the spot
            idx = floor(obj.spot_price / ds);
            idx = max(1, min(idx, Ns - 2));  % Safe indexing
            w = (obj.spot_price - sGrid(idx + 1)) / ds;
            v0 = (1 - w) * V_curr(idx) + w * V_curr(idx + 1);
        end

        function v0 = fdm_single_barrier(obj, Ns, Nt, theta, ratio, m)
            % Finite difference method for single barrier option.
            % Using a non-uniform grid, with shape controlled by input ratio.
            % Discrete monitoring.
            % Iteration process is only suitable for knock out option.
            % For knock in, this function uses a vanilla option to minus the identical knock out.
            %
            % Parameters:
            % Ns: Number of points in stock price axis
            % Nt: Number of points in time axis, per monitoring interval
            % theta:
            %     0 : Fully implicit method
            %     0.5 : Crank-nicolson method
            % ratio: The parameter used to control the shape of the grid, controls 
            % how many extra points to cluster near key values (spot and barrier)
            % 
            % m : monitoring times
            
            % Discretize Nt-1 points between every two monitoring time
            step = Nt;
            Nt = Nt * m;
            
            % Set up parameters
            mu = obj.r - obj.q;
            range = 5 * obj.sig * sqrt(obj.maturity);
            
            if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                Smax = obj.spot_price * exp((mu - obj.sig^2/2.0) * obj.maturity + range);
                Smin = obj.barrier * 0.99999999;
            elseif strcmp(obj.option_type, 'UP-AND-OUT-BARRIER') || strcmp(obj.option_type, 'UP-AND-IN-BARRIER')
                Smax = max(obj.barrier, obj.strike) * 1.0000001;
                Smin = 0;
            end
            
            % Time step
            dt = obj.maturity / Nt;
            
            % Generate non-uniform grid
            s = linspace(Smin, Smax, round(Ns * (1 - ratio) + 1));
            
            if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER')
                temp = [obj.barrier, obj.spot_price];
            elseif strcmp(obj.option_type, 'UP-AND-OUT-BARRIER')
                temp = [obj.spot_price, obj.barrier];
            elseif strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                temp = [obj.barrier, obj.spot_price];
            else
                temp = [obj.spot_price, obj.barrier];
            end
            
            % Find indices closest to key points
            lower_index = zeros(1, 2);
            upper_index = zeros(1, 2);
            for i = 1:2
                [~, idx] = min(abs(s - temp(i)));
                if s(idx) < temp(i) && idx < length(s)
                    lower_index(i) = idx;
                    upper_index(i) = idx + 1;
                else
                    lower_index(i) = max(1, idx - 1);
                    upper_index(i) = idx;
                end
            end
            
            delta_s = s(upper_index(1)) - s(lower_index(1));
            
            % Determine number of points to insert around key values
            if lower_index(2) > lower_index(1)
                count = round(Ns * ratio / 2.0);
            else
                count = round(Ns * ratio);
            end
            ds = delta_s / (count - 1);
            
            % Insert additional points around key values
            insert_vector = cell(1, 2);
            for j = 1:2
                insert_vector{j} = linspace(s(lower_index(j)) + ds, s(upper_index(j)) - ds, count);
            end
            
            % Combine grid points
            s_temp = [s(1:lower_index(1)), insert_vector{1}];
            s_temp = [s_temp, s(upper_index(1):lower_index(2))];
            
            if lower_index(2) > lower_index(1)
                s_temp = [s_temp, insert_vector{2}];
            end
            
            s_temp = [s_temp, s(upper_index(2):end)];
            s = s_temp;
            Ns = length(s) - 1;
            
            % Initialize the payoff
            if strcmp(obj.position, 'CALL')
                if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                    V_Nt = max(s - obj.strike, 0);
                    for i = 1:length(V_Nt)
                        if s(i) < obj.barrier
                            V_Nt(i) = 0;
                        end
                    end
                else
                    V_Nt = max(s - obj.strike, 0);
                    for i = 1:length(V_Nt)
                        if s(i) > obj.barrier
                            V_Nt(i) = 0;
                        end
                    end
                end
                payoff = max(s - obj.strike, 0);
            else
                if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                    V_Nt = max(obj.strike - s, 0);
                    for i = 1:length(V_Nt)
                        if s(i) < obj.barrier
                            V_Nt(i) = 0;
                        end
                    end
                else
                    V_Nt = max(obj.strike - s, 0);
                    for i = 1:length(V_Nt)
                        if s(i) > obj.barrier
                            V_Nt(i) = 0;
                        end
                    end
                end
                payoff = max(obj.strike - s, 0);
            end
            
            % Initialize the Dirichlet boundary condition
            if strcmp(obj.position, 'CALL')
                f_0 = zeros(Nt + 1, 1);
                if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                    f_Ns = Smax * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)') - ...
                           obj.strike * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)');
                else
                    f_Ns = zeros(Nt + 1, 1);
                end
            elseif strcmp(obj.position, 'PUT')
                if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                    f_0 = zeros(Nt + 1, 1);
                else
                    f_0 = obj.strike * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)');
                end
                f_Ns = zeros(Nt + 1, 1);
            end
            
            % Finite Difference Matrix Construction

            % initialize the tridiagonal matrix
            delta_s_i = 0.5 * (s(3:end) - s(1:Ns-1));
            delta_s_plus = s(3:end) - s(2:Ns);
            delta_s_minus = s(2:Ns) - s(1:Ns-1);
            
            % Set up coefficients for the matrix

            % a, b, c: Lower, diagonal, upper coefficients for implicit method.
            % alpha, beta, gamma: Same, but for the left-hand matrix.

            a = -(1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) + ...
                (1 - theta) * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            b = 1.0 / dt + (1 - theta) * obj.r + ...
                (1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) + ...
                (1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus);
            
            c = -(1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus) - ...
                (1 - theta) * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            alpha = theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) - ...
                    theta * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            beta = 1.0 / dt - theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) - ...
                   obj.r * theta - theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus);
            
            gamma = theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus) + ...
                    theta * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            % From Nt to 1, calculate V_Nt-1, V_Nt-2, ..., v0 (vectors)
            V_Nplus = V_Nt(2:Ns);
            
            % Start from maturity and go backward to time 0.
            for k = Nt:-1:1
                % V_Nplus : right-hand side of Ax=b
                V_Nplus = obj.util.dot_product(alpha, beta, gamma, V_Nplus);
                V_Nplus(1) = V_Nplus(1) - a(1) * f_0(k) + alpha(1) * f_0(k+1);
                V_Nplus(Ns-1) = V_Nplus(Ns-1) - c(Ns-1) * f_Ns(k) + gamma(Ns-1) * f_Ns(k+1);
                
                % Solve the tridiagonal system
                ab = obj.util.compress_tridiag(a, b, c);
                % Use the internally defined thomas_algorithm, not the standalone one
                V_N = obj.util.thomas_algorithm(ab, V_Nplus);
                
                % American option process
                if strcmp(obj.exercise_type, 'AMERICAN')
                    V_N = obj.projected_sor(a(2:end), b, c(1:end-1), V_Nplus, V_N, payoff(2:Ns), k, step, s(2:Ns));
                end
                
                V_Nplus = V_N;
                
                % Monitoring process
                if mod(k, step) == 0
                    if strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                        for i = 1:length(V_Nplus)
                            if s(i+1) < obj.barrier
                                V_Nplus(i) = 0;
                            end
                        end
                    else
                        for i = 1:length(V_Nplus)
                            if s(i+1) > obj.barrier
                                V_Nplus(i) = 0;
                            end
                        end
                    end
                end
            end
            
            % Linear interpolation to get price at spot price
            [~, index] = min(abs(s - obj.spot_price));
            if s(index) > obj.spot_price && index > 1
                index = index - 1;
            end
            w = (obj.spot_price - s(index)) / (s(index+1) - s(index));
            v0 = V_Nplus(index) * (1 - w) + w * V_Nplus(index+1);
            
            % Adjust for knock-in options
            if strcmp(obj.option_type, 'UP-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER')
                return;
            else
                if strcmp(obj.position, 'CALL')
                    v0 = obj.bs_call() - v0;
                else
                    if strcmp(obj.exercise_type, 'EUROPEAN')
                        v0 = obj.bs_put() - v0;
                    else
                        v0 = obj.btm_vanilla(1200) - v0;
                    end
                end
            end
        end
        
        function v0 = fdm_double_barrier(obj, Ns, Nt, theta, ratio, m)
            % Finite difference method for double barrier option.
            % Using a non-uniform grid, with shape controlled by input ratio.
            % Discrete monitoring.
            % Iteration process is only suitable for double knock out option.
            % For double knock in, this function uses a vanilla option to minus the identical knock out.
            %
            % Parameters:
            % Ns: Number of points in price axis
            % Nt: Number of points in time axis
            % theta:
            %     0 : Fully implicit method
            %     0.5 : Crank-nicolson method
            % ratio: The parameter used to control the shape of the grid
            % m : monitoring times
            
            % Discretize Nt-1 points between every two monitoring time
            step = Nt;
            Nt = Nt * m;
            
            % Set up parameters

            % mu: Drift under risk-neutral measure.
            mu = obj.r - obj.q;

            % range: Used to define the dynamic grid size (5σ width).
            range = 5 * obj.sig * sqrt(obj.maturity);
            
            % Smax, Smin: The max/min price range for the grid, 
            % slightly extended to avoid boundary issues.
            Smax = max(obj.upper_barrier, max(obj.spot_price, obj.strike) * ...
                exp((mu - obj.sig^2/2.0) * obj.maturity + range)) * 1.0000001;
            Smin = obj.lower_barrier * 0.9999999;
            
            % Time step

            % Each time step spans a fraction of the total maturity.
            dt = obj.maturity / Nt;
            
            % Generate non-uniform grid

            % creates a uniform grid for part of the domain, depending on ratio.
            s = linspace(Smin, Smax, round(Ns * (1 - ratio) + 1));
            
            % these will later receive more grid points for higher accuracy.
            temp = [obj.lower_barrier, obj.spot_price, obj.upper_barrier];
            
            % Find indices closest to key points
            lower_index = zeros(1, 3);
            upper_index = zeros(1, 3);
            for i = 1:3
                [~, idx] = min(abs(s - temp(i)));
                if s(idx) < temp(i) && idx < length(s)
                    lower_index(i) = idx;
                    upper_index(i) = idx + 1;
                else
                    lower_index(i) = max(1, idx - 1);
                    upper_index(i) = idx;
                end
            end
            
            delta_s = s(upper_index(1)) - s(lower_index(1));
            
            % Determine number of points to insert around key values
            if lower_index(2) > lower_index(1) && lower_index(3) > lower_index(2)
                count = round(Ns * ratio / 3.0);
            else
                count = round(Ns * ratio / 2.0);
            end
            ds = delta_s / (count - 1);
            
            % Insert additional points around key values
            insert_vector = cell(1, 3);
            for j = 1:3
                insert_vector{j} = linspace(s(lower_index(j)) + ds, s(upper_index(j)) - ds, count);
            end
            
            % Combine grid points
            s_temp = [s(1:lower_index(1)), insert_vector{1}];
            s_temp = [s_temp, s(upper_index(1):lower_index(2))];
            
            if lower_index(2) > lower_index(1)
                s_temp = [s_temp, insert_vector{2}];
            end
            
            s_temp = [s_temp, s(upper_index(2):lower_index(3))];
            
            if lower_index(3) > lower_index(2)
                s_temp = [s_temp, insert_vector{3}];
            end
            
            % Rebuild the s grid with the newly added concentrated points.
            s_temp = [s_temp, s(upper_index(3):end)];
            s = s_temp;
            Ns = length(s) - 1;
            
            % Initialize the payoff

            % ensures values outside the barriers are zero (knocked out).
            if strcmp(obj.position, 'CALL')
                V_Nt = max(s - obj.strike, 0);
                for i = 1:length(V_Nt)
                    if s(i) > obj.upper_barrier || s(i) < obj.lower_barrier
                        V_Nt(i) = 0;
                    end
                end
                payoff = max(s - obj.strike, 0);
            else
                V_Nt = max(obj.strike - s, 0);
                for i = 1:length(V_Nt)
                    if s(i) > obj.upper_barrier || s(i) < obj.lower_barrier
                        V_Nt(i) = 0;
                    end
                end
                payoff = max(obj.strike - s, 0);
            end
            
            % Initialize the Dirichlet boundary condition
            % Boundary conditions are zero at both ends (Dirichlet BCs), 
            % appropriate for knock-out options
            if strcmp(obj.position, "CALL") || strcmp(obj.position, "PUT")
                f_0 = zeros(Nt + 1, 1);
                f_Ns = zeros(Nt + 1, 1);
            end
            
            % Initialize the tridiagonal matrix

            % handles irregular spacing by computing finite differences 
            % for each grid cell.
            delta_s_i = 0.5 * (s(3:end) - s(1:Ns-1));
            delta_s_plus = s(3:end) - s(2:Ns);
            delta_s_minus = s(2:Ns) - s(1:Ns-1);
            
            % Set up coefficients for the matrix (time stepping)

            % They reflect the terms from the PDE (diffusion and convection) on a non-uniform grid.
            a = -(1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) + ...
                (1 - theta) * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            b = 1.0 / dt + (1 - theta) * obj.r + ...
                (1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) + ...
                (1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus);
            
            c = -(1.0 - theta) * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus) - ...
                (1 - theta) * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            alpha = theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) - ...
                    theta * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            beta = 1.0 / dt - theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_minus) - ...
                   obj.r * theta - theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus);
            
            gamma = theta * obj.sig^2 * s(2:Ns).^2 ./ (2.0 * delta_s_i .* delta_s_plus) + ...
                    theta * mu * s(2:Ns) ./ (2 * delta_s_i);
            
            % From Nt to 1, calculate V_Nt-1, V_Nt-2, ..., v0 (vectors)
            V_Nplus = V_Nt(2:Ns);
            
            for k = Nt:-1:1
                % V_Nplus : right-hand side of Ax=b
                V_Nplus = obj.util.dot_product(alpha, beta, gamma, V_Nplus);
                V_Nplus(1) = V_Nplus(1) - a(1) * f_0(k) + alpha(1) * f_0(k+1);
                V_Nplus(Ns-1) = V_Nplus(Ns-1) - c(Ns-1) * f_Ns(k) + gamma(Ns-1) * f_Ns(k+1);
                
                % Solve the tridiagonal system
                ab = obj.util.compress_tridiag(a, b, c);
                V_N = obj.util.thomas_algorithm(ab, V_Nplus);
                
                % American option process
                if strcmp(obj.exercise_type, 'AMERICAN')
                    V_N = obj.projected_sor(a(2:end), b, c(1:end-1), V_Nplus, V_N, payoff(2:Ns), k, step, s(2:Ns));
                end
                
                V_Nplus = V_N;
                
                % Monitoring process

                % every step iterations (i.e., each monitoring date), 
                % knock out values that breach barriers.
                if mod(k, step) == 0
                    for i = 1:length(V_Nplus)
                        if s(i+1) > obj.upper_barrier || s(i+1) < obj.lower_barrier
                            V_Nplus(i) = 0;
                        end
                    end
                end
            end
            
            % Linear interpolation to get price at spot price
            [~, index] = min(abs(s - obj.spot_price));
            if s(index) > obj.spot_price && index > 1
                index = index - 1;
            end
            w = (obj.spot_price - s(index)) / (s(index+1) - s(index));
            v0 = V_Nplus(index) * (1 - w) + w * V_Nplus(index+1);
            
            % Adjust for knock-in options

            % if it's a knock-out, return the price.
            if strcmp(obj.option_type, 'KNOCK-OUT-DOUBLE-BARRIER')
                return;
            else
                % if knock-in, price = vanilla - knock-out 
                % (since vanilla = knock-in + knock-out due to parity).
                if strcmp(obj.position, 'CALL')
                    v0 = obj.bs_call() - v0;
                else
                    if strcmp(obj.exercise_type, 'EUROPEAN')
                        v0 = obj.bs_put() - v0;
                    else
                        v0 = obj.btm_vanilla(1200) - v0;
                    end
                end
            end
        end
          
        function price = btm_vanilla(obj, Nt)
            % Binomial tree pricing model for Vanilla American options
            %
            % Parameters:
            % Nt: Number of time periods used in BTM model
            
            dt = obj.maturity / Nt;

            % Up and down factors
            % derived from the Cox-Ross-Rubinstein (CRR) model.

            u = exp(obj.sig * sqrt(dt));
            d = 1 / u;
            % Up-going risk-neutral probability
            p = (exp((obj.r - obj.q) * dt) - d) / (u - d);
            
            % Maturity pay-off
            s_maturity = zeros(Nt + 1, 1);
            for i = 1:Nt+1
                s_maturity(i) = obj.spot_price * u^(2*(i-1) - Nt);
            end
            
            if strcmp(obj.position, 'CALL')
                V = max(0, s_maturity - obj.strike);
            else
                V = max(0, obj.strike - s_maturity);
            end
            
            % Backward iteration
            for t = Nt-1:-1:0
                % at each node, option price = discounted expected value of future payoffs.
                V_temp = zeros(length(V)-1, 1);
                for i = 1:length(V_temp)
                    V_temp(i) = exp(-obj.r * dt) * (p * V(i+1) + (1-p) * V(i));
                end
                
                if strcmp(obj.exercise_type, 'AMERICAN')
                    % Early exercise process
                    % for CALL: compare early exercise value (S - K) with continuation value
                    if strcmp(obj.position, 'CALL')
                        for i = 1:length(V_temp)
                            V_temp(i) = max(obj.spot_price * u^(2*(i-1)-t) - obj.strike, V_temp(i));
                        end
                    else
                        % for PUT: compare early exercise value (K - S) with continuation value.
                        for i = 1:length(V_temp)
                            V_temp(i) = max(obj.strike - obj.spot_price * u^(2*(i-1)-t), V_temp(i));
                        end
                    end
                end
                % Prepare for next backward step.
                V = V_temp;
            end
            
            % final result: option price at the root of the tree (time 0, spot = S₀).
            price = V(1);
        end
        
        function x_k_plus = projected_sor(obj, array_a, array_b, array_c, b, x_0, payoff, t, step, s)
            % Projected SOR method. Used for pricing American Options under Finite Difference Framework
            %
            % Parameters:
            % array_a : lower diagonal
            % array_b : main diagonal
            % array_c : upper diagonal
            % b : right-hand side vector
            % x_0 : initial guess for the solution
            % payoff : payoff vector
            % t : current time step
            % step : monitoring step
            % s : price grid
            
            % Constructs the tridiagonal matrix A from three vectors:
            n = length(array_b);
            A = diag(array_a, -1) + diag(array_b, 0) + diag(array_c, 1);
            
            % convergence tolerance, scaled by system size.
            epsilon = 0.01 * n;

            % Relaxation factor
            % ω < 1 → under-relaxation (slow but stable)
            % ω > 1 → over-relaxation (faster convergence but risky)
            omega = 0.5;

            % x_k: the current guess for solution (starts with x_0).
            x_k = x_0;
            
            % Begin iteration: prepare new estimate x_k_plus
            while true
                x_k_plus = zeros(n, 1);
                
                % Loop over each row (grid point)
                for i = 1:n
                    % sum_1: previous value (already updated in this iteration).
                    if i == 1
                        sum_1 = 0;
                    else
                        sum_1 = A(i, i-1) * x_k_plus(i-1);
                    end
                    
                    % sum_2: next value (not yet updated, so use x_k from previous iteration).
                    if i < n
                        sum_2 = A(i, i+1) * x_k(i+1);
                    else
                        sum_2 = 0;
                    end
                    
                    % Core Gauss-Seidel update formula (isolates unknown x(i) from equation Ax = b).
                    x_k_plus_gs_temp = (-sum_1 - sum_2 + b(i)) / A(i, i);
                    x_k_plus(i) = max((1 - omega) * x_k(i) + omega * x_k_plus_gs_temp, payoff(i));
                end
                
                % Apply barrier conditions

                % if stock crosses either barrier, option is deactivated (knock-out) → value = 0.
                if strcmp(obj.option_type, 'KNOCK-OUT-DOUBLE-BARRIER') || strcmp(obj.option_type, 'KNOCK-IN-DOUBLE-BARRIER')
                    % mod(t, step) == 0 ensures checking is done only at discrete barrier monitoring points.
                    if mod(t, step) == 0
                        for i = 1:length(x_k_plus)
                            if s(i) > obj.upper_barrier || s(i) < obj.lower_barrier
                                x_k_plus(i) = 0;
                            end
                        end
                    end
                elseif strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                    if mod(t, step) == 0
                        for i = 1:length(x_k_plus)
                            if s(i) < obj.barrier
                                x_k_plus(i) = 0;
                            end
                        end
                    end
                elseif strcmp(obj.option_type, 'UP-AND-OUT-BARRIER') || strcmp(obj.option_type, 'UP-AND-IN-BARRIER')
                    if mod(t, step) == 0
                        for i = 1:length(x_k_plus)
                            if s(i) > obj.barrier
                                x_k_plus(i) = 0;
                            end
                        end
                    end
                end
                
                % Convergence check

                % check if solution has converged (using L2 norm of change).
                if sqrt(sum((abs(x_k_plus - x_k)).^2)) < epsilon
                    break;
                else
                    x_k = x_k_plus;
                end
            end
        end
   
    end
end

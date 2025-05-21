classdef option < handle
    % OPTION Class for pricing options under Black-Scholes model
    %
    % All option pricing procedures are under Black-Scholes model, with the basic assumptions:
    %   - constant risk-free rate
    %   - constant volatility
    %   - constant and continuous dividend
    
    properties
        r           % risk-free rate
        q           % continuous dividend
        spot_price  % current price of underlying
        strike      % strike price
        sig         % volatility
        maturity    % time to maturity (in years)
        option_type % type of option (VANILLA, BARRIER, etc.)
        exercise_type % exercise type (EUROPEAN, AMERICAN)
        position    % call or put
        barrier     % barrier for single barrier options
        lower_barrier % lower barrier for double barrier options
        upper_barrier % upper barrier for double barrier options
    end
    
    methods (Static)
        function p = norm_cdf(x)
            p = 0.5 * (1 + erf(x ./ sqrt(2)));
        end
    end

    methods
        function obj = option(r, q, spot_price, strike, sig, T, option_type, exercise_type, position, T_switch, varargin)
            % Constructor for option class
            %
            % Parameters:
            % r : risk free rate
            % q : continuous dividend
            % spot_price: current price of underlying
            % strike: strike price
            % sig: volatility
            % T: time to maturity
            % T_switch:
            %     Default/0 : T is counted in trading days, AKA 252 days a year.
            %     1 : T is counted in natural year days, AKA 365 days a year.
            %     2 : T is counted in year.
            % option_type : Type of the option, e.g. Vanilla, Asian, Barrier, Lookback
            % exercise_type: Exercise type of the option, e.g. European, American
            % position: call or put
            
            % Set attributes
            obj.r = r;
            obj.q = q;
            obj.spot_price = spot_price;
            obj.strike = strike;
            obj.sig = sig;
            
            % Calculate maturity in years based on T_switch
            if nargin < 10 || isempty(T_switch)
                T_switch = 0;
            end
            
            if T_switch == 0
                obj.maturity = T / 252;
            elseif T_switch == 1
                obj.maturity = T / 365;
            elseif T_switch == 2
                obj.maturity = T;
            else
                error('Invalid input of T_switch');
            end
            
            % Set option type
            valid_OptionType = {'VANILLA', 'ASIAN', 'LOOKBACK', 'UP-AND-IN-BARRIER', ...
                'UP-AND-OUT-BARRIER', 'DOWN-AND-IN-BARRIER', 'DOWN-AND-OUT-BARRIER', ...
                'KNOCK-OUT-DOUBLE-BARRIER', 'KNOCK-IN-DOUBLE-BARRIER'};
            
            option_type_upper = upper(option_type);
            if any(strcmp(option_type_upper, valid_OptionType))
                obj.option_type = option_type_upper;
            else
                error('Currently don''t support this type of options!');
            end
            
            % Set exercise type
            exercise_type_upper = upper(exercise_type);
            if strcmp(exercise_type_upper, 'EUROPEAN')
                obj.exercise_type = 'EUROPEAN';
            elseif strcmp(exercise_type_upper, 'AMERICAN')
                obj.exercise_type = 'AMERICAN';
            else
                error('Invalid Exercise Type');
            end
            
            % Set position
            position_upper = upper(position);
            if strcmp(position_upper, 'CALL')
                obj.position = 'CALL';
            elseif strcmp(position_upper, 'PUT')
                obj.position = 'PUT';
            else
                error('Invalid position type');
            end
            
            % Parse additional arguments as name-value pairs
            p = inputParser;
            p.addParameter('barrier', [], @isnumeric);
            p.addParameter('lower_barrier', [], @isnumeric);
            p.addParameter('upper_barrier', [], @isnumeric);
            p.parse(varargin{:});
            
            % Set barriers for single barrier options
            if strcmp(obj.option_type, 'UP-AND-OUT-BARRIER')
                if ~isempty(p.Results.barrier)
                    if p.Results.barrier <= obj.spot_price
                        error('Barrier should not be smaller than spot price');
                    else
                        obj.barrier = p.Results.barrier;
                    end
                else
                    error('Lack information of barrier value');
                end
            elseif strcmp(obj.option_type, 'UP-AND-IN-BARRIER')
                if ~isempty(p.Results.barrier)
                    if p.Results.barrier <= obj.spot_price
                        error('Barrier should not be smaller than spot price');
                    else
                        obj.barrier = p.Results.barrier;
                    end
                else
                    error('Lack information of barrier value');
                end
            elseif strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER')
                if ~isempty(p.Results.barrier)
                    if p.Results.barrier >= obj.spot_price
                        error('Barrier should not be larger than spot price');
                    else
                        obj.barrier = p.Results.barrier;
                    end
                else
                    error('Lack information of barrier value');
                end
            elseif strcmp(obj.option_type, 'DOWN-AND-IN-BARRIER')
                if ~isempty(p.Results.barrier)
                    if p.Results.barrier >= obj.spot_price
                        error('Barrier should not be larger than spot price');
                    else
                        obj.barrier = p.Results.barrier;
                    end
                else
                    error('Lack information of barrier value');
                end
            end
            
            % Set barriers for double barrier options
            if strcmp(obj.option_type, 'KNOCK-OUT-DOUBLE-BARRIER')
                if ~isempty(p.Results.lower_barrier) && ~isempty(p.Results.upper_barrier)
                    if p.Results.upper_barrier <= p.Results.lower_barrier
                        error('Upper barrier should be larger than lower barrier');
                    elseif obj.spot_price >= p.Results.upper_barrier || obj.spot_price <= p.Results.lower_barrier
                        error('Invalid barrier value');
                    else
                        obj.lower_barrier = p.Results.lower_barrier;
                        obj.upper_barrier = p.Results.upper_barrier;
                    end
                else
                    error('Lack information of upper/lower barrier value');
                end
            elseif strcmp(obj.option_type, 'KNOCK-IN-DOUBLE-BARRIER')
                if ~isempty(p.Results.lower_barrier) && ~isempty(p.Results.upper_barrier)
                    if p.Results.upper_barrier <= p.Results.lower_barrier
                        error('Upper barrier should be larger than lower barrier');
                    elseif obj.spot_price >= p.Results.upper_barrier || obj.spot_price <= p.Results.lower_barrier
                        error('Invalid barrier value');
                    else
                        obj.lower_barrier = p.Results.lower_barrier;
                        obj.upper_barrier = p.Results.upper_barrier;
                    end
                else
                    error('Lack information of upper/lower barrier value');
                end
            end
        end
        
        function call_price = Black_Scholes_Call(obj)
            d1 = (log(obj.spot_price / obj.strike) + (obj.r - obj.q + 0.5 * obj.sig^2) * obj.maturity) / ...
                 (obj.sig * sqrt(obj.maturity));
            d2 = d1 - obj.sig * sqrt(obj.maturity);
            forward = obj.spot_price * exp((obj.r - obj.q) * obj.maturity);
            call_price = exp(-obj.r * obj.maturity) * ...
                         (forward * option.norm_cdf(d1) - obj.strike * option.norm_cdf(d2));  % <-- call using class name
        end
    
        function put_price = Black_Scholes_Put(obj)
            d1 = (log(obj.spot_price / obj.strike) + (obj.r - obj.q + 0.5 * obj.sig^2) * obj.maturity) / ...
                 (obj.sig * sqrt(obj.maturity));
            d2 = d1 - obj.sig * sqrt(obj.maturity);
            forward = obj.spot_price * exp((obj.r - obj.q) * obj.maturity);
            put_price = exp(-obj.r * obj.maturity) * ...
                        (obj.strike * option.norm_cdf(-d2) - forward * option.norm_cdf(-d1));  % <-- call using class name
        end

        function price = BTM_Vanilla(obj, Nt)
            % Binomial tree pricing model for Vanilla American options
            %
            % Parameters:
            % Nt: Number of time periods used in BTM model
            
            dt = obj.maturity / Nt;
            % Up and down factors
            u = exp(obj.sig * sqrt(dt));
            d = 1 / u;
            % Up-going probability
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
                V_temp = zeros(length(V)-1, 1);
                for i = 1:length(V_temp)
                    V_temp(i) = exp(-obj.r * dt) * (p * V(i+1) + (1-p) * V(i));
                end
                
                if strcmp(obj.exercise_type, 'AMERICAN')
                    % Early exercise process
                    if strcmp(obj.position, 'CALL')
                        for i = 1:length(V_temp)
                            V_temp(i) = max(obj.spot_price * u^(2*(i-1)-t) - obj.strike, V_temp(i));
                        end
                    else
                        for i = 1:length(V_temp)
                            V_temp(i) = max(obj.strike - obj.spot_price * u^(2*(i-1)-t), V_temp(i));
                        end
                    end
                end
                
                V = V_temp;
            end
            
            price = V(1);
        end
        
        function ab = tri_bound_to_ab(~, coeff_m1_arr, coeff_0_arr, coeff_p1_arr)
            % Compute ab which stores the matrix elements of the tri-diagonal matrix bounded matrix
            % ab = [
            %     [0, c1, ..., c_I_3, c_I_2],
            %     [b1, b2, ...,b_I_3, b_I_1],
            %     [a2, a3, ...,a_I_1, 0]
            % ]
            
            length_arr = length(coeff_m1_arr);
            ab = zeros(3, length_arr);
            ab(1, 2:length_arr) = coeff_p1_arr(1:length_arr-1);
            ab(2, :) = coeff_0_arr(:);
            ab(3, 1:length_arr-1) = coeff_m1_arr(2:end);
        end
        
        function array_result = my_dot_product(~, array_a, array_b, array_c, vector)
            % Compute (sub-diagonal, main, super-diagonal) · vector
            %
            % array_a, array_b, array_c each length n
            % vector               length n
            %
            % result(i) = a(i)*vector(i-1) + b(i)*vector(i) + c(i)*vector(i+1)
            % with the obvious boundary checks.
    
            n = numel(vector);
            array_result = zeros(n,1);
    
            for i = 1:n
                % sub-diagonal contribution
                if i > 1
                    array_result(i) = array_result(i) + array_a(i) * vector(i-1);
                end
                % main diagonal
                array_result(i) = array_result(i) + array_b(i) * vector(i);
                % super-diagonal
                if i < n
                    array_result(i) = array_result(i) + array_c(i) * vector(i+1);
                end
            end
        end
        
        function x_k_plus = Projected_SOR(obj, array_a, array_b, array_c, b, x_0, payoff, t, step, s)
            % Projected SOR method. Used for pricing American Options under Finite Difference Framework
            %
            % Parameters:
            % array_a : -1 diagonal
            % array_b : diagonal
            % array_c : 1 diagonal
            % b : right-hand side vector
            % x_0 : initial guess
            % payoff : payoff vector
            % t : current time step
            % step : monitoring step
            % s : price grid
            
            length_arr = length(array_b);
            A = diag(array_a, -1) + diag(array_b, 0) + diag(array_c, 1);
            
            epsilon = 0.01 * length_arr;
            omega = 0.5;
            x_k = x_0;
            
            while true
                x_k_plus = zeros(length_arr, 1);
                
                for i = 1:length_arr
                    if i == 1
                        sum_1 = 0;
                    else
                        sum_1 = A(i, i-1) * x_k_plus(i-1);
                    end
                    
                    if i < length_arr
                        sum_2 = A(i, i+1) * x_k(i+1);
                    else
                        sum_2 = 0;
                    end
                    
                    x_k_plus_gs_temp = (-sum_1 - sum_2 + b(i)) / A(i, i);
                    x_k_plus(i) = max((1 - omega) * x_k(i) + omega * x_k_plus_gs_temp, payoff(i));
                end
                
                % Apply barrier conditions
                if strcmp(obj.option_type, 'KNOCK-OUT-DOUBLE-BARRIER') || strcmp(obj.option_type, 'KNOCK-IN-DOUBLE-BARRIER')
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
                if sqrt(sum((abs(x_k_plus - x_k)).^2)) < epsilon
                    break;
                else
                    x_k = x_k_plus;
                end
            end
        end
        
        function v_0 = FDM_DoubleBarrier_NonUnifromGrid(obj, Ns, Nt, theta, ratio, m)
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

            % Smax, Smin: The max/min price range for the grid, slightly extended to avoid boundary issues.
            
            range = 5 * obj.sig * sqrt(obj.maturity);
            Smax = max(obj.upper_barrier, max(obj.spot_price, obj.strike) * ...
                exp((mu - obj.sig^2/2.0) * obj.maturity + range)) * 1.0000001;
            Smin = obj.lower_barrier * 0.9999999;
            
            % Time step
            dt = obj.maturity / Nt;
            
            % Generate non-uniform grid
            s = linspace(Smin, Smax, round(Ns * (1 - ratio) + 1));
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
            
            s_temp = [s_temp, s(upper_index(3):end)];
            s = s_temp;
            Ns = length(s) - 1;
            
            % Initialize the payoff
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
            if strcmp(obj.position, "CALL")
                f_0 = zeros(Nt + 1, 1);
                f_Ns = zeros(Nt + 1, 1);
            elseif strcmp(obj.position, "PUT")
                f_0 = zeros(Nt + 1, 1);
                f_Ns = zeros(Nt + 1, 1);
            end
            
            % Initialize the tridiagonal matrix
            delta_s_i = 0.5 * (s(3:end) - s(1:Ns-1));
            delta_s_plus = s(3:end) - s(2:Ns);
            delta_s_minus = s(2:Ns) - s(1:Ns-1);
            
            % Set up coefficients for the matrix
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
            
            % From Nt to 1, calculate V_Nt-1, V_Nt-2, ..., V_0 (vectors)
            V_Nplus = V_Nt(2:Ns);
            
            for k = Nt:-1:1
                % V_Nplus : right-hand side of Ax=b
                V_Nplus = obj.my_dot_product(alpha, beta, gamma, V_Nplus);
                V_Nplus(1) = V_Nplus(1) - a(1) * f_0(k) + alpha(1) * f_0(k+1);
                V_Nplus(Ns-1) = V_Nplus(Ns-1) - c(Ns-1) * f_Ns(k) + gamma(Ns-1) * f_Ns(k+1);
                
                % Solve the tridiagonal system
                ab = obj.tri_bound_to_ab(a, b, c);
                V_N = my_thomas_algorithm(ab, V_Nplus);
                
                % American option process
                if strcmp(obj.exercise_type, 'AMERICAN')
                    V_N = obj.Projected_SOR(a(2:end), b, c(1:end-1), V_Nplus, V_N, payoff(2:Ns), k, step, s(2:Ns));
                end
                
                V_Nplus = V_N;
                
                % Monitoring process
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
            v_0 = V_Nplus(index) * (1 - w) + w * V_Nplus(index+1);
            
            % Adjust for knock-in options
            if strcmp(obj.option_type, 'KNOCK-OUT-DOUBLE-BARRIER')
                return;
            else
                if strcmp(obj.position, 'CALL')
                    v_0 = obj.Black_Scholes_Call() - v_0;
                else
                    if strcmp(obj.exercise_type, 'EUROPEAN')
                        v_0 = obj.Black_Scholes_Put() - v_0;
                    else
                        v_0 = obj.BTM_Vanilla(1200) - v_0;
                    end
                end
            end
        end
            
        function v_0 = FDM_SingleBarrier_NonUnifromGrid(obj, Ns, Nt, theta, ratio, m)
            % Finite difference method for single barrier option.
            % Using a non-uniform grid, with shape controlled by input ratio.
            % Discrete monitoring.
            % Iteration process is only suitable for knock out option.
            % For knock in, this function uses a vanilla option to minus the identical knock out.
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
            
            % Initialize the tridiagonal matrix
            delta_s_i = 0.5 * (s(3:end) - s(1:Ns-1));
            delta_s_plus = s(3:end) - s(2:Ns);
            delta_s_minus = s(2:Ns) - s(1:Ns-1);
            
            % Set up coefficients for the matrix
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
            
            % From Nt to 1, calculate V_Nt-1, V_Nt-2, ..., V_0 (vectors)
            V_Nplus = V_Nt(2:Ns);
            
            for k = Nt:-1:1
                % V_Nplus : right-hand side of Ax=b
                V_Nplus = obj.my_dot_product(alpha, beta, gamma, V_Nplus);
                V_Nplus(1) = V_Nplus(1) - a(1) * f_0(k) + alpha(1) * f_0(k+1);
                V_Nplus(Ns-1) = V_Nplus(Ns-1) - c(Ns-1) * f_Ns(k) + gamma(Ns-1) * f_Ns(k+1);
                
                % Solve the tridiagonal system
                ab = obj.tri_bound_to_ab(a, b, c);
                % Use the internally defined thomas_algorithm, not the standalone one
                V_N = my_thomas_algorithm(ab, V_Nplus);
                
                % American option process
                if strcmp(obj.exercise_type, 'AMERICAN')
                    V_N = obj.Projected_SOR(a(2:end), b, c(1:end-1), V_Nplus, V_N, payoff(2:Ns), k, step, s(2:Ns));
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
            v_0 = V_Nplus(index) * (1 - w) + w * V_Nplus(index+1);
            
            % Adjust for knock-in options
            if strcmp(obj.option_type, 'UP-AND-OUT-BARRIER') || strcmp(obj.option_type, 'DOWN-AND-OUT-BARRIER')
                return;
            else
                if strcmp(obj.position, 'CALL')
                    v_0 = obj.Black_Scholes_Call() - v_0;
                else
                    if strcmp(obj.exercise_type, 'EUROPEAN')
                        v_0 = obj.Black_Scholes_Put() - v_0;
                    else
                        v_0 = obj.BTM_Vanilla(1200) - v_0;
                    end
                end
            end
        end
        
        function v_0 = FDM_Vanilla_Implicit(obj, Ns, Nt, m)
            % Finite difference method for vanilla option.
            % Trivial implicit method.
            %
            % Parameters:
            % Ns: Number of points in price axis
            % Nt: Number of points in time axis
            % m : monitoring times
            
            % Discretize Nt-1 points between every two monitoring time
            step = Nt;
            Nt = Nt * m;
            
            % Set up parameters
            mu = obj.r - obj.q;
            range = 3 * obj.sig * sqrt(obj.maturity);
            Smax = obj.spot_price * exp((mu - obj.sig^2/2.0) * obj.maturity + range);
            Smin = 0;
            
            % Time step
            dt = obj.maturity / Nt;
            ds = (Smax - Smin) / Ns;  % totally Ns + 1 in column grid
            
            % Grid points
            sGrid = linspace(Smin, Smax, Ns + 1);
            
            % Initialize the payoff
            if strcmp(upper(obj.position), 'CALL')
                V_Nt = max(sGrid - obj.strike, 0);
            elseif strcmp(upper(obj.position), 'PUT')
                V_Nt = max(obj.strike - sGrid, 0);
            end
            
            % Initialize payoff for all grid points
            if strcmp(upper(obj.position), 'CALL')
                payoff = max(sGrid - obj.strike, 0);
            else
                payoff = max(obj.strike - sGrid, 0);
            end
            
            % Initialize the Dirichlet boundary condition
            if strcmp(upper(obj.position), 'CALL')
                f_0 = zeros(Nt + 1, 1);
                f_Ns = Smax * exp(-obj.q * linspace(0, obj.maturity, Nt + 1)') - ...
                       obj.strike * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)');
            elseif strcmp(upper(obj.position), 'PUT')
                f_0 = obj.strike * exp(-obj.r * linspace(0, obj.maturity, Nt + 1)');
                f_Ns = zeros(Nt + 1, 1);
            else
                error('Invalid option_type!!');
            end
            
            % Initialize the tridiagonal matrix by scalar-form
            i = (1:Ns-1)';
            
            % From a_2 to a_I-1 are in the calculation matrix
            a = -(obj.sig^2 * i.^2 - (obj.r - obj.q) * i) * dt / 2.0;
            
            % From b_1 to b_I-1 are in the calculation matrix
            b = 1 + obj.sig^2 * i.^2 * dt + obj.r * dt;
            
            % From c_1 to c_I-2 are in the calculation matrix
            c = -(obj.sig^2 * i.^2 + (obj.r - obj.q) * i) * dt / 2.0;
            
            % From Nt to 1, calculate V_Nt-1, V_Nt-2, ..., V_0 (vectors)
            V_Nplus = V_Nt(2:Ns);
            
            for k = Nt:-1:1
                % Update boundary conditions for the current time step
                V_Nplus(1) = V_Nplus(1) - a(1) * f_0(k);
                V_Nplus(Ns-1) = V_Nplus(Ns-1) - c(Ns-1) * f_Ns(k);
                
                % Construct the tridiagonal matrix
                ab = obj.tri_bound_to_ab(a, b, c);
                % Use the internally defined thomas_algorithm, not the standalone one
                V_N = my_thomas_algorithm(ab, V_Nplus);
                
                % American option process
                if strcmp(obj.exercise_type, 'AMERICAN')
                    V_N = obj.Projected_SOR(a(2:end), b, c(1:end-1), V_Nplus, V_N, payoff(2:Ns), k, step, sGrid(2:Ns));
                end
                
                V_Nplus = V_N;
            end
            
            % Linear interpolation to get price at spot price
            index = floor(obj.spot_price / ds);
            index = max(1, min(index, Ns-1)); % Ensure index is within bounds
            
            w = (obj.spot_price - sGrid(index+1)) / (sGrid(index+2) - sGrid(index+1));
            w = max(0, min(w, 1)); % Ensure weight is between 0 and 1
            
            v_0 = V_N(index) * (1 - w) + w * V_N(index+1);
        end
        
        function price = Monte_Carlo_Vanilla(obj, path_num)
            % Monte Carlo method for European vanilla option.
            %
            % Parameter:
            % path_num : number of simulation times
            
            mu = obj.r - obj.q;
            simulated_prices = zeros(path_num, 1);
            
            for i = 1:path_num
                simulated_prices(i) = obj.spot_price * exp((mu - 0.5 * obj.sig^2) * obj.maturity + ...
                                      obj.sig * sqrt(obj.maturity) * randn);
            end
            
            if strcmp(obj.position, 'CALL')
                simulated_option_prices = max(simulated_prices - obj.strike, 0);
            else
                simulated_option_prices = max(obj.strike - simulated_prices, 0);
            end
            
            price = mean(simulated_option_prices) * exp(-obj.r * obj.maturity);
        end
    end
end

% The thomas_algorithm function needs to be defined as a nested function
% within the option class methods or as a local function within the file.
% Here's a corrected implementation:
function output = my_thomas_algorithm(ab, d)
    % Thomas algorithm for tridiagonal systems
    % ab: [upper diagonal; main diagonal; lower diagonal]
    % d: right hand side vector
    
    n = length(d);
    c_prime = zeros(n, 1);
    d_prime = zeros(n, 1);
    
    % Forward sweep
    c_prime(1) = ab(1, 2) / ab(2, 1);
    d_prime(1) = d(1) / ab(2, 1);
    
    for i = 2:n-1
        c_prime(i) = ab(1, i+1) / (ab(2, i) - ab(3, i) * c_prime(i-1));
    end
    
    for i = 2:n
        d_prime(i) = (d(i) - ab(3, i) * d_prime(i-1)) / (ab(2, i) - ab(3, i) * c_prime(i-1));
    end
    
    % Back substitution
    x = zeros(n, 1);
    x(n) = d_prime(n);
    
    for i = n-1:-1:1
        x(i) = d_prime(i) - c_prime(i) * x(i+1);
    end
    
    output = x;
end
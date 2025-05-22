classdef AnalyticalSolution
    % AnalyticalSolution - Black-Scholes closed-form solution for European options

    properties
        sigma       % Volatility
        r           % Risk-free interest rate
        K           % Strike price
        T           % Time to maturity
        optionType  % 'call' or 'put'
    end

    methods (Static)
        function p = norm_cdf(x)
            % Computes the cumulative distribution function of standard normal distribution
            p = 0.5 * (1 + erf(x ./ sqrt(2)));
        end
    end

    methods
        function obj = AnalyticalSolution(sigma, r, K, T, optionType)
            obj.sigma = sigma;
            obj.r = r;
            obj.K = K;
            obj.T = T;
            obj.optionType = lower(optionType);  % Ensure lowercase
        end

        function price = solve(obj, S, t)
            % Returns price at time t (i.e., time to maturity is T - t)
            tau = obj.T - t;  % Time to maturity

            if tau <= 0
                % Option is at maturity
                if strcmp(obj.optionType, 'call')
                    price = max(S - obj.K, 0);
                else
                    price = max(obj.K - S, 0);
                end
                return;
            end

            % Standard Black-Scholes formula
            d1 = (log(S ./ obj.K) + (obj.r + 0.5 * obj.sigma^2) * tau) / (obj.sigma * sqrt(tau));
            d2 = d1 - obj.sigma * sqrt(tau);

            if strcmp(obj.optionType, 'call')
                price = S .* AnalyticalSolution.norm_cdf(d1) - ...
                        obj.K * exp(-obj.r * tau) .* AnalyticalSolution.norm_cdf(d2);
            else
                price = obj.K * exp(-obj.r * tau) .* AnalyticalSolution.norm_cdf(-d2) - ...
                        S .* AnalyticalSolution.norm_cdf(-d1);
            end
        end
    end
end

classdef utils
    methods (Static)
        
        function p = norm_cdf(x)
            % Computes the cumulative distribution function of standard normal distribution
            p = 0.5 * (1 + erf(x ./ sqrt(2)));
        end

        function ab = compress_tridiag(a, b, c)
            % Constructs a compact 3xN representation of a tridiagonal matrix
            %
            % a: lower diagonal (a_2 to a_n)
            % b: main diagonal (b_1 to b_n)
            % c: upper diagonal (c_1 to c_{n-1})
            % ab: 3xN matrix representation [upper; main; lower]

            n = length(b);
            ab = zeros(3, n);
            ab(1, 2:n) = c(1:n-1);           % upper diagonal
            ab(2, :)    = b(:);          % main diagonal
            ab(3, 1:n-1) = a(2:end);     % lower diagonal (shifted)
        end

        function res = dot_product(a, b, c, v)
            % Efficient dot product for tridiagonal matrix T * v
            %
            % a, b, c: diagonals (length n)
            % v: vector (length n)
            % res: result vector T*v
            
            n = numel(v);
            res = zeros(n, 1);

            res(1)     = b(1)*v(1) + c(1)*v(2);
            for i = 2:n-1
                res(i) = a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1);
            end
            res(n)     = a(n)*v(n-1) + b(n)*v(n);
        end

        function x = thomas_algorithm(A, d)
            % Solves a tridiagonal system Ax = d using Thomas algorithm
            %
            % A: 3xN compressed matrix [upper; main; lower]
            % d: right-hand side vector
            % x: solution vector
            
            n = length(d);
            c = zeros(n, 1); % modified upper diag
            d_ = zeros(n, 1); % modified rhs

            % Forward sweep
            c(1) = A(1,2) / A(2,1);
            d_(1) = d(1) / A(2,1);

            for i = 2:n-1
                denom = A(2,i) - A(3,i) * c(i-1);
                c(i) = A(1,i+1) / denom;
                d_(i) = (d(i) - A(3,i) * d_(i-1)) / denom;
            end
            d_(n) = (d(n) - A(3,n) * d_(n-1)) / (A(2,n) - A(3,n) * c(n-1));

            % Back substitution
            x = zeros(n,1);
            x(n) = d_(n);
            for i = n-1:-1:1
                x(i) = d_(i) - c(i) * x(i+1);
            end
        end

    end
end

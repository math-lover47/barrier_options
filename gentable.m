function gentable(strikes, analytical_vals, fdm_vals, option_type, timer_start)
    % TABLE Generates a formatted output table comparing analytical and FDM values
    %
    % Inputs:
    %   strikes          - Vector of strike prices
    %   analytical_vals  - Vector of analytical option values
    %   fdm_vals         - Vector of FDM-computed option values
    %   option_type      - 'call' or 'put'
    %   timer_start      - time value from tic (to compute total time)

    if nargin < 5
        timer_start = tic;  % start timer if not provided
    end

    % Header
    if strcmpi(option_type, 'call')
        fprintf('\nTest for call options\n\n');
    elseif strcmpi(option_type, 'put')
        fprintf('\nTest for put options\n\n');
    else
        error('Invalid option_type. Use "call" or "put".');
    end

    % Table Header
    fprintf('%-14s %-20s %-20s\n', ...
        'strike', 'analytical', 'fdm');

    % Table Body
    for i = 1:length(strikes)
        fprintf('%-14.0f %-20.10f %-20.10f\n', ...
            strikes(i), analytical_vals(i), fdm_vals(i));
    end

    % Time Output
    elapsed = toc(timer_start);
    fprintf('\n------------over----------\n');
    fprintf('Total time taken: %.4f seconds\n\n', elapsed);
end

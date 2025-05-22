function run()
    % Main runner function for Black-Scholes European Option pricing
    % Compares Analytical, FEM P1, and FEM P2 solutions
    clc; close all;
    fprintf('=== Black-Scholes European Option Pricing ===\n');
    
    % Run test case for call option
    solveBSEquation('optionType', 'call', 'K', 100, 'T', 0.25, ...
                   'sigma', 0.2, 'r', 0.05, 'nx', 200, 'nt', 50, ...
                   'S_max', 200, 'visualize', true, 'plot3D', true);

    fprintf('\n=== Analysis completed ===\n');
end

function solutions = solveBSEquation(varargin)
    % Main solver function with Analytical, FEM P1, and FEM P2 solutions
    try
        params = parseInputs(varargin);
        [xvec, tvec, analytical, fem_p1, fem_p2] = setupDomain(params);
        
        % Compute solutions
        analytical_sol = computeAnalyticalSolution(analytical, xvec, tvec);
        fem_p1_sol = computeFEMSolution(fem_p1, params, tvec, xvec);
        fem_p2_sol = computeFEMSolution(fem_p2, params, tvec, xvec);
        
        % Package and visualize
        solutions = packageSolutions(analytical_sol, fem_p1_sol, fem_p2_sol, xvec, tvec, params);
        if params.visualize
            visualizeSolutions2D(xvec, tvec, solutions);
            if params.plot3D
                visualizeSolutions3D(xvec, tvec, solutions);
            end
        end
        
        % Compute and display errors
        computeAndDisplayErrors(solutions);
        
    catch ME
        fprintf('Error: %s\n', ME.message);
        rethrow(ME);
    end
end

function [xvec, tvec, analytical, fem_p1, fem_p2] = setupDomain(params)
    % Initialize analytical and FEM solvers
    xvec = linspace(0, params.S_max, params.nx)';
    tvec = linspace(0, params.T, params.nt);
    
    % Analytical solution
    analytical = AnalyticalSolution(params.sigma, params.r, params.K, params.T, params.optionType);
    
    % FEM P1 solver
    fem_p1 = FEMSolver(params.nx-1, 'P1', params.S_max, params.theta, ...
                      params.K, params.sigma, params.r, params.T, params.optionType);
    
    % FEM P2 solver (use fewer elements since P2 has more DOF per element)
    fem_p2 = FEMSolver(floor((params.nx-1)/2), 'P2', params.S_max, params.theta, ...
                      params.K, params.sigma, params.r, params.T, params.optionType);
end

function analytical_sol = computeAnalyticalSolution(analytical, xvec, tvec)
    % Compute analytical solution
    analytical_sol = zeros(length(xvec), length(tvec));
    
    for i = 1:length(tvec)
        analytical_sol(:, i) = analytical.solve(xvec, tvec(i));
    end
end

function fem_sol = computeFEMSolution(fem, params, tvec, xvec)
    % Compute and interpolate FEM solution
    fem_sol_raw = zeros(fem.dof, length(tvec));
    fem.set_payoff_initial_condition();
    fem_sol_raw(:, end) = fem.u0;
    
    dt = params.T / (length(tvec)-1);
    for i = length(tvec)-1:-1:1
        fem.step(dt);
        fem_sol_raw(:, i) = fem.u0;
    end
    
    % Interpolate to common grid
    fem_sol = zeros(length(xvec), length(tvec));
    for i = 1:length(tvec)
        fem_sol(:, i) = interp1(fem.nodes, fem_sol_raw(:, i), xvec, 'spline', 'extrap');
    end
end

function visualizeSolutions2D(xvec, tvec, solutions)
    % Create 2D comparison plots
    figure('Position', [100 100 1400 900]);
    
    % Current time comparison (t = 0)
    subplot(2,3,1);
    plot(xvec, solutions.analytical(:,1), 'k-', 'LineWidth', 2.5); hold on;
    plot(xvec, solutions.fem_p1(:,1), 'b--', 'LineWidth', 2);
    plot(xvec, solutions.fem_p2(:,1), 'r:', 'LineWidth', 2);
    title('Option Value at t = 0 (Current Time)', 'FontSize', 12);
    legend({'Analytical', 'FEM P1', 'FEM P2'}, 'Location', 'best', 'FontSize', 10);
    xlabel('Asset Price (S)', 'FontSize', 10);
    ylabel('Option Value', 'FontSize', 10);
    grid on;
    
    % Mid-time comparison
    mid_idx = round(length(tvec)/2);
    subplot(2,3,2);
    plot(xvec, solutions.analytical(:,mid_idx), 'k-', 'LineWidth', 2.5); hold on;
    plot(xvec, solutions.fem_p1(:,mid_idx), 'b--', 'LineWidth', 2);
    plot(xvec, solutions.fem_p2(:,mid_idx), 'r:', 'LineWidth', 2);
    title(sprintf('Option Value at t = %.3f', tvec(mid_idx)), 'FontSize', 12);
    legend({'Analytical', 'FEM P1', 'FEM P2'}, 'Location', 'best', 'FontSize', 10);
    xlabel('Asset Price (S)', 'FontSize', 10);
    ylabel('Option Value', 'FontSize', 10);
    grid on;
    
    % Payoff at maturity (t = T)
    subplot(2,3,3);
    plot(xvec, solutions.analytical(:,end), 'k-', 'LineWidth', 2.5); hold on;
    plot(xvec, solutions.fem_p1(:,end), 'b--', 'LineWidth', 2);
    plot(xvec, solutions.fem_p2(:,end), 'r:', 'LineWidth', 2);
    title('Payoff at Maturity (t = T)', 'FontSize', 12);
    legend({'Analytical', 'FEM P1', 'FEM P2'}, 'Location', 'best', 'FontSize', 10);
    xlabel('Asset Price (S)', 'FontSize', 10);
    ylabel('Payoff', 'FontSize', 10);
    grid on;
    
    % Error plots
    subplot(2,3,4);
    error_p1 = abs(solutions.fem_p1(:,1) - solutions.analytical(:,1));
    error_p2 = abs(solutions.fem_p2(:,1) - solutions.analytical(:,1));
    semilogy(xvec, error_p1, 'b--', 'LineWidth', 2); hold on;
    semilogy(xvec, error_p2, 'r:', 'LineWidth', 2);
    title('Absolute Error at t = 0', 'FontSize', 12);
    legend({'FEM P1 Error', 'FEM P2 Error'}, 'Location', 'best', 'FontSize', 10);
    xlabel('Asset Price (S)', 'FontSize', 10);
    ylabel('Absolute Error', 'FontSize', 10);
    grid on;
    
    % Time evolution for specific asset price
    S_strike_idx = find(xvec >= solutions.params.K, 1);
    subplot(2,3,5);
    plot(tvec, solutions.analytical(S_strike_idx,:), 'k-', 'LineWidth', 2.5); hold on;
    plot(tvec, solutions.fem_p1(S_strike_idx,:), 'b--', 'LineWidth', 2);
    plot(tvec, solutions.fem_p2(S_strike_idx,:), 'r:', 'LineWidth', 2);
    title(sprintf('Time Evolution at S = %.1f (≈ Strike)', xvec(S_strike_idx)), 'FontSize', 12);
    legend({'Analytical', 'FEM P1', 'FEM P2'}, 'Location', 'best', 'FontSize', 10);
    xlabel('Time (t)', 'FontSize', 10);
    ylabel('Option Value', 'FontSize', 10);
    grid on;
    
    % Price evolution comparison (multiple times)
    subplot(2,3,6);
    time_indices = round(linspace(1, length(tvec), 4));
    colors = lines(4);
    for i = 1:length(time_indices)
        idx = time_indices(i);
        plot(xvec, solutions.analytical(:,idx), '-', 'Color', colors(i,:), 'LineWidth', 2); hold on;
        plot(xvec, solutions.fem_p1(:,idx), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
        plot(xvec, solutions.fem_p2(:,idx), ':', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    title('Price Evolution Over Time', 'FontSize', 12);
    xlabel('Asset Price (S)', 'FontSize', 10);
    ylabel('Option Value', 'FontSize', 10);
    grid on;
    legend_str = arrayfun(@(i) sprintf('t=%.3f', tvec(i)), time_indices, 'UniformOutput', false);
    legend(legend_str, 'Location', 'best', 'FontSize', 9);
    
    sgtitle('Black-Scholes Solution Comparison', 'FontSize', 14, 'FontWeight', 'bold');
end

function visualizeSolutions3D(xvec, tvec, solutions)
    % Create 3D surface plots
    figure('Position', [200 200 1400 500]);
    [T, S] = meshgrid(tvec, xvec);
    
    % Analytical solution
    subplot(1,3,1);
    surf(S, T, solutions.analytical, 'EdgeColor', 'none');
    title('Analytical Solution', 'FontSize', 12);
    xlabel('Asset Price (S)'); ylabel('Time (t)'); zlabel('Option Value');
    view(45, 30);
    colorbar;
    
    % FEM P1 solution
    subplot(1,3,2);
    surf(S, T, solutions.fem_p1, 'EdgeColor', 'none');
    title('FEM P1 Solution', 'FontSize', 12);
    xlabel('Asset Price (S)'); ylabel('Time (t)'); zlabel('Option Value');
    view(45, 30);
    colorbar;
    
    % FEM P2 solution
    subplot(1,3,3);
    surf(S, T, solutions.fem_p2, 'EdgeColor', 'none');
    title('FEM P2 Solution', 'FontSize', 12);
    xlabel('Asset Price (S)'); ylabel('Time (t)'); zlabel('Option Value');
    view(45, 30);
    colorbar;
    
    sgtitle('3D Solution Surfaces', 'FontSize', 14, 'FontWeight', 'bold');
end

function computeAndDisplayErrors(solutions)
    % Compute and display error metrics
    fprintf('\n=== Error Analysis ===\n');
    
    % Compute errors at current time (t = 0)
    error_p1_current = abs(solutions.fem_p1(:,1) - solutions.analytical(:,1));
    error_p2_current = abs(solutions.fem_p2(:,1) - solutions.analytical(:,1));
    
    % Compute RMS errors
    rms_p1 = sqrt(mean(error_p1_current.^2));
    rms_p2 = sqrt(mean(error_p2_current.^2));
    
    % Compute max errors
    max_p1 = max(error_p1_current);
    max_p2 = max(error_p2_current);
    
    % Display results
    fprintf('At Current Time (t = 0):\n');
    fprintf('  FEM P1 - RMS Error: %.6e, Max Error: %.6e\n', rms_p1, max_p1);
    fprintf('  FEM P2 - RMS Error: %.6e, Max Error: %.6e\n', rms_p2, max_p2);
    
    % Relative errors at strike price
    S_strike_idx = find(solutions.xvec >= solutions.params.K, 1);
    analytical_at_strike = solutions.analytical(S_strike_idx, 1);
    
    if analytical_at_strike > 1e-10
        rel_error_p1 = abs(solutions.fem_p1(S_strike_idx,1) - analytical_at_strike) / analytical_at_strike;
        rel_error_p2 = abs(solutions.fem_p2(S_strike_idx,1) - analytical_at_strike) / analytical_at_strike;
        
        fprintf('At Strike Price (S = %.1f):\n', solutions.xvec(S_strike_idx));
        fprintf('  Analytical Value: %.6f\n', analytical_at_strike);
        fprintf('  FEM P1 Value: %.6f (Rel. Error: %.4f%%)\n', ...
                solutions.fem_p1(S_strike_idx,1), rel_error_p1*100);
        fprintf('  FEM P2 Value: %.6f (Rel. Error: %.4f%%)\n', ...
                solutions.fem_p2(S_strike_idx,1), rel_error_p2*100);
    end
    
    fprintf('======================\n');
end

function solutions = packageSolutions(analytical_sol, fem_p1_sol, fem_p2_sol, xvec, tvec, params)
    % Package solutions into structured output
    solutions = struct();
    solutions.analytical = analytical_sol;
    solutions.fem_p1 = fem_p1_sol;
    solutions.fem_p2 = fem_p2_sol;
    solutions.xvec = xvec;
    solutions.tvec = tvec;
    solutions.params = params;
    solutions.timestamp = datetime('now');
end

function params = parseInputs(args)
    % Parse input arguments with default values
    p = inputParser;
    
    % Add parameters with validation
    addParameter(p, 'optionType', 'call', @(x) ismember(lower(x), {'call', 'put'}));
    addParameter(p, 'K', 100, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'T', 0.25, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'sigma', 0.2, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'r', 0.05, @(x) isnumeric(x) && x >= 0);
    addParameter(p, 'nx', 100, @(x) isnumeric(x) && x > 10);
    addParameter(p, 'nt', 50, @(x) isnumeric(x) && x > 10);
    addParameter(p, 'S_max', 200, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'theta', 0.5, @(x) isnumeric(x) && x >= 0 && x <= 1);
    addParameter(p, 'visualize', true, @islogical);
    addParameter(p, 'plot3D', false, @islogical);
    
    parse(p, args{:});
    params = p.Results;
    params.optionType = lower(params.optionType);
    
    displayParameters(params);
end

function displayParameters(params)
    % Display parameters in formatted way
    fprintf('\n--- Parameters ---\n');
    fprintf('  Option Type: %s\n', upper(params.optionType));
    fprintf('  Strike Price (K): %.2f\n', params.K);
    fprintf('  Time to Maturity (T): %.4f\n', params.T);
    fprintf('  Volatility (σ): %.4f\n', params.sigma);
    fprintf('  Risk-free Rate (r): %.4f\n', params.r);
    fprintf('  Grid Points: %d × %d (S × t)\n', params.nx, params.nt);
    fprintf('  Max Asset Price: %.2f\n', params.S_max);
    fprintf('  Theta (time scheme): %.2f\n', params.theta);
    fprintf('------------------\n');
end
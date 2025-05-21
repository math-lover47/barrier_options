function visualization(S, T, V_analytical, V_fdm)
    % VISUALIZATION Plots analytical and FDM results for option pricing
    % 
    % Inputs:
    %   S - Stock price grid (vector, length Ns+1)
    %   T - Time grid (vector, length Nt+1)
    %   V_analytical - Matrix of analytical values (size Nt+1 x Ns+1)
    %   V_fdm        - Matrix of FDM solution (size Nt+1 x Ns+1)
    
    [S_grid, T_grid] = meshgrid(S, T);

    % 1. Plot Analytical Surface
    figure;
    subplot(2,2,1);
    surf(S_grid, T_grid, V_analytical, 'EdgeColor', 'none');
    title('Analytical Solution (Black-Scholes)');
    xlabel('Stock Price');
    ylabel('Time to Maturity');
    zlabel('Option Value');
    colormap jet;
    view(45, 30);
    colorbar;
    grid on;

    % 2. Plot FDM Surface
    subplot(2,2,2);
    surf(S_grid, T_grid, V_fdm, 'EdgeColor', 'none');
    title('Finite Difference Solution');
    xlabel('Stock Price');
    ylabel('Time to Maturity');
    zlabel('Option Value');
    colormap turbo;
    view(45, 30);
    colorbar;
    grid on;

    % 3. Compare Solutions at Time = 0 (i.e., T(end))
    subplot(2,1,2);
    plot(S, V_analytical(end,:), 'b-', 'LineWidth', 2); hold on;
    plot(S, V_fdm(end,:), 'r--', 'LineWidth', 2);
    title('Comparison at Time Zero');
    xlabel('Stock Price');
    ylabel('Option Value');
    legend('Analytical', 'FDM', 'Location', 'Best');
    grid on;

end

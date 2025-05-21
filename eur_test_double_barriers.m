clc;
clear;
close all;

% --- SHARED PARAMETERS ---
r = 0.05;
q = 0.0;
sigma = 0.2;
T_years = 10;
strike = 100;
pos = 'CALL';                 % 'CALL' or 'PUT'
exercise = 'EUROPEAN';        % Only EUROPEAN supported for now
theta = 0.5;                  % Crank-Nicolson
t_switch = 2;                 % Not used here, placeholder
m = 1;                        % Monitoring frequency

% --- DOUBLE BARRIER CONFIG ---
barrier_type = 'KNOCK-IN-DOUBLE-BARRIER';  % <- change to KNOCK-IN-... to test knock-in
upper_barrier = 130;
lower_barrier = 70;

% --- GRID SETUP ---
S = linspace(lower_barrier + 1, upper_barrier - 1, 100);   % Spot price range inside barriers
T = linspace(0, T_years, 100);                             % Time points for visualization

% --- PRICING PARAMETERS ---
Ns = 80;       % Grid points for price
Nt = 100;      % Grid points for time

% --- ALLOCATE STORAGE ---
V_fdm = zeros(length(T), length(S));
V_bs = zeros(size(V_fdm));     % For comparison (when applicable)

% --- LOOP OVER SPOT PRICES ---
for i = 1:length(S)
    spot = S(i);

    % Create double barrier option
    opt_barrier = option_new(r, q, spot, strike, sigma, T_years, ...
                             barrier_type, exercise, pos, ...
                             t_switch, theta, ...
                             'lower_barrier', lower_barrier, ...
                             'upper_barrier', upper_barrier);

    % Create vanilla option (for knock-in reference)
    opt_vanilla = option_new(r, q, spot, strike, sigma, T_years, ...
                             'VANILLA', exercise, pos, ...
                             t_switch, theta);

    % --- PRICE WITH FDM ---
    Vbar = opt_barrier.fdm_double_barrier(Ns, Nt, theta, 0.3, m);
    Vvan = 0;

    % --- KNOCK-IN: Use vanilla - knock-out ---
    if contains(barrier_type, 'IN')
        if strcmp(pos, 'CALL')
            Vvan = opt_vanilla.bs_call();
        else
            Vvan = opt_vanilla.bs_put();
        end
        V_bs(end, i) = Vvan - Vbar;
    else
        % Knock-OUT value is directly computed
        V_bs(end, i) = NaN;   % No BS comparison here
    end

    V_fdm(end, i) = Vbar;
end

% --- VISUALIZATION ---
if contains(barrier_type, 'IN')
    % Show knock-in comparison
    visualization(S, T, repmat(V_bs(end, :), length(T), 1), repmat(V_fdm(end, :), length(T), 1));
else
    % Only FDM (no BS model for knock-out)
    visualization(S, T, NaN(size(V_fdm)), repmat(V_fdm(end, :), length(T), 1));
end

% --- PRINT TABLE ---
timer_start = tic;
gentable(S, V_bs(end, :), V_fdm(end, :), lower(pos), timer_start);

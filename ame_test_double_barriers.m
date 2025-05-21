clc;
clear;
close all;

% --- SHARED PARAMETERS ---
r = 0.05;
q = 0.0;
sigma = 0.65;
T_years = 0.5;
strike = 0.5;
pos = 'PUT';                  % Choose 'CALL' or 'PUT'
exercise = 'AMERICAN';        % American option test
theta = 0.5;                  % Crank-Nicolson scheme
t_switch = 2;                 % Placeholder (not used here)
m = 1;                        % Monitoring frequency (1 = annual)

% --- DOUBLE BARRIER CONFIGURATION ---
barrier_type = 'KNOCK-OUT-DOUBLE-BARRIER';   % Knock-IN or Knock-OUT
lower_barrier = 0.4;                         % Must be < spot
upper_barrier = 0.65;                         % Must be > spot

% --- GRID SETUP ---
S = linspace(lower_barrier + 0.00001, upper_barrier - 0.00001, 100);  % Stock prices
T = linspace(0, T_years, 100);                                % Time axis for visualization

% --- FDM PARAMETERS ---
Ns = 100;         % Space grid
Nt = 100;        % Time grid per monitoring
ratio = 0.2;     % Concentration ratio for non-uniform grid

% --- ALLOCATE STORAGE ---
V_fdm = zeros(length(T), length(S));
V_bs = zeros(size(V_fdm));    % For knock-in: vanilla - knockout

% --- MAIN LOOP OVER SPOT PRICES ---
for i = 1:length(S)
    spot = S(i);

    % Barrier option object
    opt_barrier = option_new(r, q, spot, strike, sigma, T_years, ...
                             barrier_type, exercise, pos, ...
                             t_switch, theta, ...
                             'lower_barrier', lower_barrier, ...
                             'upper_barrier', upper_barrier);
    
    % Vanilla option object for knock-in handling
    opt_vanilla = option_new(r, q, spot, strike, sigma, T_years, ...
                             'VANILLA', exercise, pos, ...
                             t_switch, theta);

    % --- Price using FDM ---
    Vbar = opt_barrier.fdm_double_barrier(Ns, Nt, theta, ratio, m);
    Vvan = 0;

    % --- Knock-in handling ---
    if contains(barrier_type, 'IN')
        if strcmp(pos, 'CALL')
            Vvan = opt_vanilla.bs_call();
        else
            if strcmp(exercise, 'EUROPEAN')
                Vvan = opt_vanilla.bs_put();
            else
                Vvan = opt_vanilla.btm_vanilla(1200);  % American vanilla using binomial tree
            end
        end
        V_bs(end, i) = Vvan - Vbar;  % Knock-in = Vanilla - Knock-out
    else
        % Knock-out directly modeled by FDM
        V_bs(end, i) = 0;  % Optional: analytical benchmark not available
    end

    V_fdm(end, i) = Vbar;
end

% --- VISUALIZATION ---
visualization(S, T, repmat(V_bs(end, :), length(T), 1), repmat(V_fdm(end, :), length(T), 1));

% --- RESULT TABLE ---
timer_start = tic;
gentable(S, V_bs(end, :), V_fdm(end, :), lower(pos), timer_start);

clc;
clear;
close all;

% --- SHARED PARAMETERS ---
r = 0.05;
q = 0.0;
sigma = 0.2;
T_years = 1;
strike = 200;
pos = 'CALL';                 % 'CALL' or 'PUT'
exercise = 'AMERICAN';        % Only change for AMERICAN test
theta = 0.5;                  % Crank-Nicolson scheme
t_switch = 2;                 % Not used in barrier pricing but kept for interface consistency
m = 1;                        % Monitoring frequency (e.g., 1 for annual, 12 for monthly)

% --- BARRIER CONFIGURATION ---
barrier_type = 'DOWN-AND-OUT-BARRIER';  % Change as needed
barrier_level = 180;                    % Must be < spot for DOWN, > spot for UP

% --- GRID SETUP ---
S = linspace(barrier_level + 0.1, 250, 100);   % Stock prices for evaluation
T = linspace(0, T_years, 100);                 % Time steps (for visualization)

% --- PRICING PARAMETERS ---
Ns = 80;        % Spatial grid resolution
Nt = 100;        % Temporal resolution per monitoring interval

% --- ALLOCATE STORAGE ---
V_fdm = zeros(length(T), length(S));
V_bs = zeros(size(V_fdm));     % Will be replaced by BTM or BS if needed

% --- MAIN LOOP OVER SPOT PRICES ---
for i = 1:length(S)
    spot = S(i);

    % Create barrier option object
    opt_barrier = option_new(r, q, spot, strike, sigma, T_years, ...
                             barrier_type, exercise, pos, ...
                             t_switch, theta, 'barrier', barrier_level);

    % Create vanilla option for knock-in if needed
    opt_vanilla = option_new(r, q, spot, strike, sigma, T_years, ...
                             'VANILLA', exercise, pos, ...
                             t_switch, theta);

    % --- Finite Difference Pricing ---
    Vbar = opt_barrier.fdm_single_barrier(Ns, Nt, theta, 0.2, m);
    Vvan = 0;

    % --- Knock-In Logic: use vanilla - knock-out
    if contains(barrier_type, 'IN')
        if strcmp(pos, 'CALL')
            Vvan = opt_vanilla.bs_call();
        else
            if strcmp(exercise, 'EUROPEAN')
                Vvan = opt_vanilla.bs_put();
            else
                Vvan = opt_vanilla.btm_vanilla(1200);  % American vanilla using binomial
            end
        end
        V_bs(end, i) = Vvan - Vbar;
    else
        % Knock-out is directly computed
        V_bs(end, i) = 0;  % Optional to skip direct analytical comparison
    end

    V_fdm(end, i) = Vbar;
end

% --- VISUALIZATION ---
visualization(S, T, repmat(V_bs(end, :), length(T), 1), repmat(V_fdm(end, :), length(T), 1));

% --- PRINT RESULT TABLE ---
timer_start = tic;
gentable(S, V_bs(end, :), V_fdm(end, :), lower(pos), timer_start);

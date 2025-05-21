clc;
clear;
close all;

% --- SHARED PARAMETERS ---
r = 0.05;
q = 0.0;
sigma = 0.2;
T_years = 1;
strike = 200;
pos = 'CALL';
exercise = 'EUROPEAN';
theta = 0.5;
t_switch = 2;   % Time is in years
m = 1;

% --- BARRIER CONFIG ---
barrier_type = 'DOWN-AND-OUT-BARRIER';  % <- change this for other types
barrier_level = 90;                     % For DOWN barriers: < spot price

% --- GRID SETUP ---
S = linspace(barrier_level + 1, 120, 200);    % Stock prices for evaluation
T = linspace(0, T_years, 100); % Time steps

% --- PRICING PARAMETERS ---
Ns = 50;      % Stock grid resolution
Nt = 100;      % Time grid resolution

% --- ALLOCATE STORAGE ---
V_fdm = zeros(length(T), length(S));
V_bs = zeros(size(V_fdm));

% --- LOOP OVER SPOT PRICES ---
for i = 1:length(S)
    spot = S(i);
    
    % Create barrier option object
    opt_barrier = option_new(r, q, spot, strike, sigma, T_years, ...
                             barrier_type, exercise, pos, ...
                             t_switch, theta, 'barrier', barrier_level);

    % Create vanilla option for comparison
    opt_vanilla = option_new(r, q, spot, strike, sigma, T_years, ...
                             'VANILLA', exercise, pos, ...
                             t_switch, theta);

    % --- PRICE WITH FDM ---
    Vbar = opt_barrier.fdm_single_barrier(Ns, Nt, theta, 0.2, m);
    Vvan = 0;

    % --- PRICE WITH BLACK-SCHOLES for Knock-IN or Knock-OUT logic ---
    if contains(barrier_type, 'IN')
        % Knock-IN = Vanilla - Knock-OUT
        if strcmp(pos, 'CALL')
            Vvan = opt_vanilla.bs_call();
        else
            Vvan = opt_vanilla.bs_put();
        end
        V_bs(end, i) = Vvan - Vbar;
    else
        % Knock-OUT is directly modeled
        V_bs(end, i) = 0;  % Optional: skip BS comparison for knock-out
    end

    V_fdm(end, i) = Vbar;
end

% --- VISUALIZATION ---
visualization(S, T, repmat(V_bs(end, :), length(T), 1), repmat(V_fdm(end, :), length(T), 1));

% --- PRINT TABLE ---
timer_start = tic;
gentable(S, V_bs(end, :), V_fdm(end, :), lower(pos), timer_start);

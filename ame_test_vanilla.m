clc;
clear;
close all;

%% --- SHARED PARAMETERS ---
r = 0.05;                 % risk-free rate
q = 0.0;                  % dividend yield
sigma = 0.2;              % volatility
T_years = 1;              % time to maturity
strike = 100;             % strike price
pos = 'PUT';              % position: CALL or PUT
type = 'VANILLA';         % option type
exercise = 'AMERICAN';    % American option
theta = 0.5;              % Crank-Nicolson (theta = 0.5)
t_switch = 2;             % dummy, not used here
m = 1;                    % sub-stepping per time step (for FDM)

%% --- STOCK PRICE AND TIME GRIDS ---
S = linspace(50, 150, 1000);   % Stock price grid
T = linspace(0, T_years, 100); % Time grid

%% --- INITIALIZE AMERICAN OPTION OBJECT ---
opt_am = option_new(r, q, 100, strike, sigma, T_years, ...
                    type, exercise, pos, t_switch, theta);

V_fdm = zeros(length(T), length(S));

%% --- GRID PARAMETERS FOR FDM ---
Ns = 10;        % Number of spatial steps
Nt = 500;        % Number of time steps

%% --- COMPUTE FDM PRICES AT T = 0 ---
for i = 1:length(S)
    opt_am.spot_price = S(i);
    V_fdm(end, i) = opt_am.fdm_vanilla(Ns, Nt, m, theta);
end

%% --- VISUALIZE FDM SOLUTION ---
visualization(S, T, repmat(V_fdm(end, :), length(T), 1), repmat(V_fdm(end, :), length(T), 1));

%% --- TABLE AT T = 0 ---
timer_start = tic;
gentable(S, V_fdm(end, :), V_fdm(end, :), 'PUT', timer_start);

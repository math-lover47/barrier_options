clc;
clear;
close all;

% --- SHARED PARAMETERS ---
r = 0.05;
q = 0.0;
sigma = 0.2;
T_years = 1;
strike = 100;
pos = 'CALL';
type = 'VANILLA';
exercise = 'EUROPEAN';
theta = 0.5;
t_switch = 2;  % years
m = 1;

% --- STOCK PRICE AND TIME GRIDS ---
S = linspace(0, 100, 1000);    % 1000 stock prices
T = linspace(0, T_years, 100); % 100 time steps

% --- EXAMPLE 1: ACCURATE FDM MATCH ---
opt1 = option_new(r, q, 100, strike, sigma, T_years, ...
                  type, exercise, pos, t_switch, theta);

V_bs = zeros(length(T), length(S));
V_fdm = zeros(size(V_bs));

Ns = 6; Nt = 1000;

for i = 1:length(S)
    opt1.spot_price = S(i);
    V_bs(end, i) = opt1.bs_call();
    V_fdm(end, i) = opt1.fdm_vanilla(Ns, Nt, m, theta);
end

% --- VISUALIZE ---
visualization(S, T, repmat(V_bs(end, :), length(T), 1), repmat(V_fdm(end, :), length(T), 1));

% --- TABLE AT T = 0 ---
timer_start = tic;
gentable(S, V_bs(end, :), V_fdm(end, :), 'call', timer_start);

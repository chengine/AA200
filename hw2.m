close all
clear all

%% a
% Parameters
global H Re Uinf
H = 2;
Re = 50000;
Uinf = 100;
cp_min = -1.5;
x_rec = 0.21;

% Solve for laminar BL
% Solve for initial conditions to turbulent layer
% cp_min happens at the beginning of the turbulent layer
% because turbulent layer has adverse pressure gradient
ue_init = sqrt((1-cp_min)*Uinf^2);
% use laminar BL to solve for theta_init
Rex = Re/Uinf * ue_init * x_rec;
theta_init = 0.671*x_rec / sqrt(Rex);
y0 = [theta_init; ue_init];

% Solve for states
N_turb = 200;
tspan = linspace(x_rec, 1, N_turb);
[T, Y] = ode45(@(t, x) deriv(t, x), tspan, y0);
theta = Y(:, 1);
ue = Y(:, 2);

figure(1)
plot(T, ue)
title('U_e(x) on upper recovery')
xlabel('X/c')

N_lam = 50;
lam_span = linspace(0, x_rec - x_rec/N_lam, N_lam);
ue_lam = ue_init*ones(N_lam, 1);
ue = [ue_lam; ue];
cp = 1 - (ue / Uinf).^2;
uel = 2*Uinf - ue;
cpl = 1 - (uel/Uinf).^2;
T = [lam_span'; T];

figure(2)
plot(T, cp)
hold on
plot(T, cpl)
hold on
xline(x_rec)
hold off
ylabel('Cp')
xlabel('X/c')
title('Camberline Pressures')

%% b and c
%Inverse design
gam = (cpl - cp)/2;
x = T;
w = 0*x;
for i = 1:length(x)
    x_temp = x;
    x_temp(i) = [];
    gam_temp = gam;
    gam_temp(i) = [];
    Y = gam_temp./ (x(i) - x_temp);
    w(i) = 1/(2*pi) * trapz(x_temp, Y);
end
alpha = trapz(x, w);
dzdx = alpha - w;   
z = cumtrapz(x, dzdx);

figure(3)
%plot((1-cos(th))/2, c)
plot(x, z)
title('Camberline')
yline(0)
xlabel('X/c')
ylabel('Z')
xline(x_rec)
axis equal

% Calculate Cl and Cm
cl = trapz(x, gam*2);
cm = trapz(x, 2*gam.*x);

%% d
Re = 5000000;
x_rec = 0.869;

% Solve for laminar BL
% Solve for initial conditions to turbulent layer
% cp_min happens at the beginning of the turbulent layer
% because turbulent layer has adverse pressure gradient
ue_init = sqrt((1-cp_min)*Uinf^2);
% use laminar BL to solve for theta_init
Rex = Re/Uinf * ue_init * x_rec;
theta_init = 0.671*x_rec / sqrt(Rex);
y0 = [theta_init; ue_init];

% Solve for states
N_turb = 1000;
tspan = linspace(x_rec, 1, N_turb);
[T, Y] = ode45(@(t, x) deriv(t, x), tspan, y0);
theta = Y(:, 1);
ue = Y(:, 2);

figure(4)
plot(T, ue)
title('U_e(x) on upper recovery')
xlabel('X/c')

N_lam = 1000;
lam_span = linspace(0, x_rec - x_rec/N_lam, N_lam);
ue_lam = ue_init*ones(N_lam, 1);
ue = [ue_lam; ue];
cp = 1 - (ue / Uinf).^2;
uel = 2*Uinf - ue;
cpl = 1 - (uel/Uinf).^2;
T = [lam_span'; T];

figure(5)
plot(T, cp)
hold on
plot(T, cpl)
hold on
xline(x_rec)
hold off
ylabel('Cp')
title('Camberline Pressures')
xlabel('X/c')

%Inverse design to solve cp, cm, and camber
gam = (cpl - cp)/2;
x = T;
w = 0*x;
for i = 1:length(x)
    x_temp = x;
    x_temp(i) = [];
    gam_temp = gam;
    gam_temp(i) = [];
    Y = gam_temp./ (x(i) - x_temp);
    w(i) = 1/(2*pi) * trapz(x_temp, Y);
end

alpha_2 = trapz(x, w);
dzdx = alpha - w;
z = cumtrapz(x, dzdx);

figure(6)
plot(x, z)
title('Camberline')
yline(0)
xlabel('X/c')
ylabel('Z')
xline(x_rec)
axis equal

% Calculate Cl and Cm
cl_2 = trapz(x, gam*2);
cm_2 = trapz(x, 2*gam.*x);
%% Functions
function dxdt = deriv(t, x)
    global H Re Uinf
    theta = x(1);
    ue = x(2);
    
    H1 = 1.5501*(H-0.6778)^(-3.064) + 3.3;
    F = 0.0306*(H1-3)^(-0.6169);
    rho_nu = Re/Uinf;
    Retheta = rho_nu*ue*theta;
    Cf = (0.246*10^(-0.678*H))/(Retheta^(0.268));
    
    dthetadt = (Cf/2 - (H+2)*(F/H1))/(1 - (H+2));
    dxdt = [dthetadt; (F/H1 - dthetadt)*ue/theta];
end

function S = create_fourier(theta)
    N = size(theta, 1);
    S = zeros(N);
    for a = 1:N
        S(a, :) = create_series(theta(a, :), N);
    end
end

function s = create_series(theta, N)
    s = cot(theta/2);
    for a = 1:N-1
        s = [s, sin(a*theta)];
    end
end

function [A, S] = find_cof(theta, cp)
    S = create_fourier(theta);
    A = inv(S)*cp;
end

function [c, alpha, th] = find_camberline(cof, Npts)
    N = size(cof, 1);
    func = @(n, theta) (cos(n*theta).*cos(theta) + n*sin(n*theta).*sin(theta))/(n^2 - 1);
    alpha = cof(1) - (1/2)*cof(2)*(-0.5*(cos(pi)^2 - cos(0)^2));
    for a = 2:N-1
        alpha = alpha - (1/2)*cof(a + 1)*(func(a, pi) - func(a, 0));
    end
    func2 = @(theta) -cos(theta);
    th = linspace(0 + pi/Npts, pi - pi/Npts, Npts);
    z = zeros(1, Npts);
    c = (1/2)*(alpha - cof(1))*(func2(th) - func2(z)) + (1/2)*cof(2)*(-0.5*(cos(th).^2 - cos(z).^2));
    for b = 2:N-1
        c = c + (1/2)*cof(b+1)*(func(b, th) - func(b, z));
    end
end
    
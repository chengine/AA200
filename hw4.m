close all
clear all

%% Q1
%Parameters
m = 35;
g = 9.81;
u0 = 15;
LD = 10;
pmax = 1000;
bmax = 12;
rho = 0.122;

N = 10000;
epsilon = 0.1;
r = linspace(epsilon, 6, N);

%For a non-accelerating aircraft, T = D, and L = G (force of gravity)
D = m*g / 10;
A = pi*r.^2;
u = (-u0 + sqrt(u0^2 + 4*D./(2*rho*A)))/2;

Pin = 2*rho*A.*(u0 + u).^2 .*u;
ind = Pin < pmax;
Pin = Pin(ind);

eta = (D*u0)./Pin;

real_r = r(ind);

figure(1)
plot(real_r, Pin)
xlabel('Radius (m)')
title('Inviscid Power')
xlim([real_r(1), real_r(end)]);

figure(2)
plot(real_r, eta)
xlabel('Radius (m)')
title('Efficiency')
xlim([real_r(1), real_r(end)]);

%% Q2a

R = real_r(1667);
p_inv = Pin(1667);
eta_inv = eta(1667);
omega = 100000;

[p,t,etaprop,cp,ct,u,lambda,r,incidence,chord,cl] = ...
    AA200OptProp2022(0.1,R,2,20,omega,u0,rho,D, 0, 1);

%% Q2b
n_iter = 10;
omega = linspace(64.7, 64.75, n_iter);
power = zeros(n_iter, 1);
parfor i=1:n_iter
    w = omega(i);
    [p,t,etaprop,cp,ct,u,lambda,r,incidence,chord,cl] = ...
        AA200OptProp2022(0.1,R,2,20,w,u0,rho,D, 1, 0);
    power(i) = p;
end
figure(5)
plot(omega, power)

%% Q2b optimal
w_opt = 64.7;
[p,t,etaprop,cp,ct,u,lambda,r,incidence,chord,cl] = ...
    AA200OptProp2022(0.1,R,2,20,w_opt,u0,rho,D, 1, 1);
N = 10000;
epsilon = 0.0000001;
cpl = [0.21, 0.21, 0.02];
cpu = [-0.33, -0.33, 0.02];
x = [0.02, 0.8, 1];

cpdiff = cpl - cpu;

% Find shape
int_gamma = @(a, b, x1, x2, y) -a*(x2-x1) - (b + a*y).*log(abs(1 - (x1-x2)./(x1-y)));
int_gamma_whole = @(y) int_gamma(0.54/0.02, 0, 0, 0.02, y) + int_gamma(0, 0.54, 0.02, 0.8, y) + int_gamma(-0.54/0.2, 0.8*(0.54/0.2)+0.54, 0.8, 1, y);
y_pts = linspace(0+epsilon, 1-epsilon, N);
z_eval_pts = int_gamma_whole(y_pts);
shape = -1/(4*pi)*cumtrapz(y_pts, z_eval_pts);
A = [y_pts(1), 1; y_pts(end), 1];
F = [-shape(1);-shape(end)];
params = inv(A)*F;
geo = params(1)*y_pts + shape + params(end);

%Find lift coefficient
alpha = params(1);
pi_pts = linspace(0+epsilon, pi-epsilon, N);
dzdx = (1/(4*pi^2))*int_gamma_whole((1 - cos(pi_pts))/2);
a0 = cumtrapz(pi_pts, dzdx);
a0 = a0(end);

dzdxcos1 = int_gamma_whole((1 - cos(pi_pts))/2).* cos(pi_pts);
a1 = -1/(2*pi^2)*cumtrapz(pi_pts, dzdxcos1);
a1 = a1(end);

dzdxcos2 = int_gamma_whole((1 - cos(pi_pts))/2).* cos(2*pi_pts);
a2 = -1/(2*pi^2)*cumtrapz(pi_pts, dzdxcos2);
a2 = a2(end);

%Lift coefficient
cl = 2*pi*(a0 + a1/2);

%Pitching moment about leading edge
cm = -cl/4 - pi/4 * (a1 - a2);

% Find center of pressure
poly = @(a, b, x1, x2) a*x2^3 / 3 + b*x2^2 / 2 - (a*x1^3 / 3 + b*x1^2 / 2);
lin = @(a, b, x1, x2) a*x2^2 / 2 + b*x2 - (a*x1^2 / 2 + b*x1);
numer = poly(0.54/0.02, 0, 0, 0.02) + poly(0, 0.54, 0.02, 0.8) + poly(-0.54/0.2, 0.8*(0.54/0.2) + 0.54, 0.8, 1);
denom = lin(0.54/0.02, 0, 0, 0.02) + lin(0, 0.54, 0.02, 0.8) + lin(-0.54/0.2, 0.8*(0.54/0.2) + 0.54, 0.8, 1);

center_press = numer/denom;

figure(1)
plot(y_pts, geo)
yline(0)
title('Camberline')
ylabel('Z')
xlabel('X in fractions of chord')
axis equal
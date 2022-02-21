close all
clear all
%% Part a

%Parameters: Q1a
aroot = 8;
atip = 8;
sweep = 0;
b = 10;
croot = 1;
ctip = 1;
Sref = 10;
u0 = 10;
rho = 1;

npanels = 100;

[yctl,gamma,chord,Cl,CL,CDi,e,Crm] = ...
  Weissinger2022(b,croot,ctip,aroot,atip,npanels,rho,u0);

figure(1)
plot(2*yctl/b, rho*u0*gamma)
title('Lift per unit span')
xlabel('\eta')

figure(2)
plot(2*yctl/b, Cl)
title('Local lift coefficient')
xlabel('\eta')

%% Part b
CL_b = CL;
CDi_b = CDi;
e_b = e;
Crm_b = Crm;
bending = 1/4 * rho * u0^2 *Sref * CL * b/4 * Crm;

%% Part c
N = 100;
epsilon = 0.001;
[same_CL, e_compare, maximum, aR_opt, aT_opt, CDi_compare, minimum, aR_optCd, aT_optCd] =...
    constArea_changeAoA(N, b,croot,ctip,npanels,rho,u0, CL_b, epsilon);

best_e = maximum;
best_Cdi = minimum;
best_alpha = [aR_opt, aT_opt];
best_alphaCd = [aR_optCd, aT_optCd];
    
%% Part 2a
%Test wing span variety
b_span = b*linspace(0.9, 1.2, 5);
b_span = [b_span(1:2), b, b_span(3:end)];
%But enforce that area is constant
c_span = 2*Sref ./ b_span;
c_span = c_span/2;

epsilon = 0.001;
Cd_store = zeros(1, 6);
e_store = zeros(1, 6);
for j = 1:6
    b_test = b_span(j);
    c_test = c_span(j);
    [same_CL, e_compare, maximum, aR_opt, aT_opt, CDi_compare, minimum, aR_optCd, aT_optCd] =...
    constArea_changeAoA(N, b_test,c_test,c_test,npanels,rho,u0, CL_b, epsilon);
    Cd_store(:, j) = minimum;
    e_store(:, j) = maximum;
end

%% Plotting 2a
prandtl = CL_b^2 ./ (pi * e_store.* b_span./c_span);
figure(3)
plot(b_span, Cd_store);
hold on
plot(b_span, prandtl)
hold off
title('Minimum induced drag with varying span');
xlabel('Span');
legend('Numerical', 'Prandtl')

figure(4)
plot(b_span, e_store);
title('Span efficiency with varying span');
xlabel('Span');
    
%% Part 2b
npanels = 200;
epsilon = 0.001;
N = 300;
Cd2_store = zeros(1, 6);
e2_store = zeros(1, 6);
for j = 1:6
    b_test = b_span(j);
    c_test = c_span(j);
    [same_CL, e_compare, maximum, aR_opt, aT_opt, CDi_compare, minimum, aR_optCd, aT_optCd] =...
    constAreaLiftBend_changeAoA(N, b_test,c_test,c_test,npanels,rho,u0, CL_b, bending, epsilon, Sref);
    Cd2_store(:, j) = minimum;
    e2_store(:, j) = maximum;
end

%% Plotting 2b
figure(5)
plot(b_span, Cd2_store)
title('Minimum induced drag with varying span');
xlabel('Span');

figure(6)
plot(b_span, e2_store);
title('Span efficiency with varying span');
xlabel('Span');

%Induced Drag
drag = 1/2 *rho*u0^2 * Sref * Cd2_store;
drag_ref = drag(3);
drag = drag / drag_ref;
figure(7)
plot(b_span, drag)
title('Minimum induced drag (normalized) with varying span');
xlabel('Span');

%% Functions
function [same_CL, e_compare, maximum, aR_opt, aT_opt, CDi_compare, minimum, aR_optCd, aT_optCd] =...
    constArea_changeAoA(N, b,croot,ctip,npanels,rho,u0, CL_compare, epsilon)
    alphaR_test = linspace(8.5, 10.5, N); 
    alphaT_test = linspace(4.5, 6.5, N);

    CL_test = zeros(N);
    CDi_test = zeros(N);
    e_test = zeros(N);

%     CL_test = zeros(1,N^2);
%     CDi_test = zeros(1,N^2);
%     e_test = zeros(1,N^2);

    for m = 1:N
        %Changing root AoA
        aR = alphaR_test(m);
        parfor k = 1:N
            %Changing tip AoA
            aT = alphaT_test(k);
            [~,~,~,~,CL,CDi,e,~] = ...
               Weissinger2022(b,croot,ctip,aR,aT,npanels,rho,u0);
           CL_test(m,k) = CL;
           CDi_test(m,k) = CDi;
           e_test(m,k) = e;
%         CL_test(:,m) = CL;
%         CDi_test(:,m) = CDi;
%         e_test(:,m) = e;
%         disp(m)
        end
        disp(m)
    end
%     CL_test = reshape(CL_test, [N, N])';
%     CDi_test = reshape(CDi_test, [N, N])';
%     e_test = reshape(e_test, [N, N])';

    % Find similar CL and find max span efficiency
    Diff = abs(CL_test - CL_compare);
    same_CL = Diff < epsilon;
    e_compare = e_test.*same_CL;
    CDi_compare = CDi_test.*same_CL;
    CDi_compare(CDi_compare == 0) = Inf;
    
    % Find maximum span efficiency
    maximum = max(max(e_compare));
    [x,y]=find(e_compare==maximum);
    aR_opt = alphaR_test(x);
    aT_opt = alphaT_test(y);
    
    % Find minimum induced drag
    minimum = min(min(CDi_compare));
    [xCd,yCd]=find(CDi_compare==minimum);
    aR_optCd = alphaR_test(xCd);
    aT_optCd = alphaT_test(yCd);
end

function [same_CL, e_compare, maximum, aR_opt, aT_opt, CDi_compare, minimum, aR_optCd, aT_optCd] =...
    constAreaLiftBend_changeAoA(N, b,croot,ctip,npanels,rho,u0, CL_compare, bending, epsilon, Sref)
    alphaR_test = linspace(0, 15, N); 
    alphaT_test = linspace(0, 15, N);

    CL_test = zeros(N);
    CDi_test = zeros(N);
    Crm_test = zeros(N);
    e_test = zeros(N);

    %Changing root AoA
    for m = 1:N
        aR = alphaR_test(m);
        %Changing tip AoA
        parfor k = 1:N
            aT = alphaT_test(k);
            [~,~,~,~,CL,CDi,e,Crm] = ...
               Weissinger2022(b,croot,ctip,aR,aT,npanels,rho,u0);
           CL_test(m,k) = CL;
           Crm_test(m,k) = Crm;
           CDi_test(m,k) = CDi;
           e_test(m,k) = e;
        end
        disp(m)
    end

    % Find similar CL and Crm and find max span efficiency
    Diff = abs(CL_test - CL_compare);
    same_CL = Diff < epsilon;
    e_compare = e_test.*same_CL;
    CDi_compare = CDi_test.*same_CL;

    Diff_crm = abs(1/4 * rho * u0^2 *Sref * CL_compare * b/4 *Crm_test - bending);
    same_Crm = Diff_crm < (1/4 * rho * u0^2 *Sref * CL_compare * b/4 * epsilon);
    e_compare = e_compare.*same_Crm;
    CDi_compare = CDi_compare.*same_Crm;
    CDi_compare(CDi_compare == 0) = Inf;
    
    % Find maximum span efficiency
    maximum = max(max(e_compare));
    if maximum == 0
        error('Unable to satisfy tolerance')
    end
    [x,y]=find(e_compare==maximum);
    aR_opt = alphaR_test(x);
    aT_opt = alphaT_test(y);
    
    % Find minimum induced drag
    minimum = min(min(CDi_compare));
    if minimum == Inf
        error('Unable to satisfy tolerance')
    end
    [xCd,yCd]=find(CDi_compare==minimum);
    aR_optCd = alphaR_test(xCd);
    aT_optCd = alphaT_test(yCd);
end
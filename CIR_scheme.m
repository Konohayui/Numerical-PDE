clc; clear all; close all;

%% Initialization
c = 335; rho = 7000/c^2;
dx = 0.005; dt = 1.5e-5;
x0 = -0.75; xf = 0.75;
T = 0.0015; t = 0:dt:T;
r = dt/dx; M = T/dt;
I = eye(2,2);

%% Initial Solution
U0 = inic(x0, xf, M, dx, c, 0);
U = U0;

%% Exact Solution
ur = inic(x0, xf, 0, dx, -c, T);
ul = inic(x0, xf, 0, dx, c, T);
g1r = ur(1,:); g1l = ul(1,:);
g2r = ur(2,:); g2l = ul(2,:);
u1_exact = (g1r + g1l + c*rho*g2r - c*rho*g2l)/2;
u2_exact = (g1r - g1l + rho*c*g2r + rho*c*g2l)/(2*c*rho);

%% Main Calculation
lambdap = zeros(2,2);
lambdam = lambdap;
jmin = x0 - M*dx; jmax = xf + M*dx;

for m = 1:M
    jmin = jmin + dx; jmax = jmax - dx;
    n = length(jmin:dx:jmax);
    new_U = zeros(2,n);
    
    for j = 2:n+1 
%         A = [U(2,j-1) U(1,j-1); c^2/U(1,j-1) U(2,j-1)];
        A = [0 c^2*rho; 1/rho 0];
        [V, lambda] = eig(A.');
        Vt = V.';
        lambdap = abs(max(lambda, 0));
        lambdam = abs(min(lambda, 0));
        d = (I - r*abs(lambda))*Vt*U(:,j) + ...
            (lambdap*Vt*U(:,j-1) + lambdam*Vt*U(:,j+1))*r;
        new_U(:,j-1) = Vt\d;
    end
    U = new_U;
    
    subplot(2,1,1)
    plot(x0:dx:xf,u1_exact,'r-');
    hold on
    plot(jmin:dx:jmax,U(1,:),'b-','markerfacecolor','b');
    hold off
    title(['Numerical and exact solution of U_1 at step ',num2str(m)])
    subplot(2,1,2)
    plot(x0:dx:xf,u2_exact,'b-')
    hold on
    plot(jmin:dx:jmax,U(2,:),'r-','markerfacecolor','r');
    hold off
    title(['Numerical and exact solution of U_2 at step ',num2str(m)])
    shg
    pause(dt);
    
end

%% Plot Figure
figure
plot(jmin:dx:jmax, U(1,:),'r-')
hold on
plot(jmin:dx:jmax,u1_exact,'b-')
hold off
legend('Numerical','Exact')
title('Numerical and Exact Solution of U_1')

figure
plot(jmin:dx:jmax, U(2,:),'r-')
hold on
plot(jmin:dx:jmax,u2_exact,'b-')
hold off
legend('Numerical','Exact')
title('Numerical and Exact Solution of U_2')

%% Initial Condition
function U0 = inic(x0, xf, M, dx, c, t)
x = x0-M*dx:dx:xf+M*dx;
l = length(x);
g1 = zeros(1,l);

for k = 1:l
    if abs(x(k) + c*t) <= 0.1
        g1(k) = 7200;
    else
        g1(k) = 7000;
    end
end

g2 = zeros(1,l);
U0 = [g1; g2];
end
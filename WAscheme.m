clear all; close all;
% function U = WAscheme(N, M, l, T, mu, D, IC, BC0, BC1)
IC =@(x) exp(x); BC0 =@(t) exp(t); BC1 =@(t) exp(1+t);
% IC =@(x) 1+1000*sin(pi*x)+1e-9*sin(5*pi*x); 
% BC0 =@(t) 1 + 0*t; 
% BC1 =@(t) 1 + 0*t;
% IC =@(x) 100 + 0*x; BC0 =@(t) 0*t; BC1 =@(t) 0*t;
% Crank Nicolson: mu =1/2
% Explicit: mu = 0
% Implicit: mu = 1

N = 30; M = 1000; l = 1; T = 2; mu = 1; D = 1;
%
% 
% initialization 
dx = l/(N+1); dt = T/M;
x = 0:dx:l; 
ini = IC(x); 
r = dt/dx^2;
fprintf(['r is ', num2str(r), '\n'])

% construct matrix A and B
Aupper = (-mu*r*D)*ones(N-1, 1);
Acenter = (1+2*mu*r*D)*ones(N, 1);
Alower = Aupper;
A = diag(Aupper, 1) + diag(Acenter) + diag(Alower, -1);
%fprintf(['Determinant of A is ', num2str(det(A)), '\n'])

Bupper = (1-mu)*r*D*ones(N-1, 1);
Bcenter = (1 -2*(1-mu)*r*D)*ones(N, 1);
Blower = Bupper; 
B = diag(Bupper, 1) + diag(Bcenter) + diag(Blower, -1);
%fprintf(['Determinant of B is ', num2str(det(B)), '\n'])

% main calculation
U = ini(2:end-1)';
for ite = 1:M
    R = B*U;
    R(1) = R(1) + r*D*mu*BC0(ite*dt) + (1-mu)*r*D*BC0((ite - 1)*dt); 
    R(end) = R(end) + r*mu*D*BC1(ite*dt) + (1 - mu)*r*D*BC1((ite - 1)*dt);
    U = A\R;
end

exact =@(x) exp(x + T);
% exact =@(x) 1+1000*sin(pi*x)*exp(-pi^2*T)+(1e-9)*sin(5*pi*x)*exp(-25*pi^2*T);
% exact =@(x) (400/pi)*sin(pi*x)*exp(-0.1*pi^2*T) + ...
%     (400/(3*pi))*sin(3*pi*x)*exp(-0.1*9*pi^2*T) + ...
%     (400/(5*pi))*sin(3*pi*x)*exp(-0.1*25*pi^2*T);
max_error = max(abs(exact(x(2:end-1))' - U))
% end

% When N = 40 and M = 500, the maximum error is 7.2201e-05 with mu = 1/2
% When N = 50 and M = 60000, the maximum error is 2.1808e-05 with mu = 0
% When N = 100 and M = 20000, the maximum error is 8.2528e-05 with mu = 1

% find approximation
X = x(2:end-1);
Y = U;
plot(X, Y)

LineH = get(gca, 'Children');
X = get(LineH, 'XData');
Y = get(LineH, 'YData');
xx = 1; 
yy = interp1(X, Y, xx);  
xx2 = X(find(abs(x - 0.6)<0.1));  
appro = 0;
L = length(xx2); 
for n = 1:L
    index = find(X == xx2(n)); 
    appro = appro + Y(index);
end
disp(appro/L)
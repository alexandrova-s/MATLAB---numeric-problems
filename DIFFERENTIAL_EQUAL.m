f = @(t, y) [y(2); -1/3 * y(2) - 1/4 * y(1)];
y0 = [0; 1];

options = odeset('RelTol', 100*eps, 'AbsTol', eps);
h = 0.01;
t = 0 : h : 10;
n = length(t);

y = radau(h, t);
[t_ode, y_ode] = ode113(f, [0,10], y0, options);

Yeuler = euler_(h, n);


figure;
plot(t, y, 'b', t_ode, y_ode(:, 1), 'k', t, Yeuler, 'r');
legend('Radau IIA', 'ode113', 'Euler');
xlabel('t');
ylabel('y');
title('Wykres przyblizenia (h = 0.01), oraz rozwiazania');


H = logspace(-4, 0.5, 20);
n = length(H);
radau2 = zeros(1, n);
radauInf = zeros(1, n);
euler2 = zeros(1, n);
eulerInf = zeros(1, n);


for k = 1 : n   
   t = 0:H(k):10;
   yRadau = radau(H(k), t);
   
   [t_ode, y_ode1] = ode113(f, t, y0, options);
   y_ode = y_ode1(:, 1);
   radau2(k) = norm(yRadau - y_ode', 2) / norm(y_ode, 2);
   radauInf(k) = norm(yRadau - y_ode', inf) / norm(y_ode, inf);
   
   yEuler = euler_(H(k), length(t));
   euler2(k) = norm(yEuler - y_ode', 2) / norm(y_ode, 2);
   eulerInf(k) = norm(yEuler - y_ode', inf) / norm(y_ode, inf);
end



figure;
loglog(H, radau2, H, euler2);
grid on;
title('Wykres bledu w normie kwadratowej');
xlabel('h');
ylabel('blad sredniokwadratowy');
legend('Radau IIA', 'Euler');

figure;
loglog(H, radauInf, H, eulerInf);
grid on;
title('Wykres bledu w normie maksymowej');
xlabel('h');
ylabel('blad maksymalny');
legend('Radau', 'Euler');


function Y = radau(h, t)

n = length(t);
Y = zeros(1, n);
Yy = zeros(2, n);
Yy(:,1) = [0; 1];
I = eye(2);

A = [ 0, 1; ...
     -1/4, -1/3];

M = [  11/45 - (7 * sqrt(6))/360,         37/225 - (169 * sqrt(6))/1800,      -2/225 + sqrt(6)/75; ...
     37/225 + (169 * sqrt(6))/1800,         11/45 + (7 * sqrt(6))/360,          -2/225 - sqrt(6)/75; ...
            4/9 - sqrt(6)/36,                    4/9 + sqrt(6)/36,                     1/9];

AA = [ I-M(1,1)*h*A,    -M(1,2)*h*A,     -M(1,3)*h*A; ...
       -M(2,1)*h*A,      I-M(2,2)*h*A,   -M(2,3)*h*A; ...
       -M(3,1)*h*A,      -M(3,2)*h*A,     I-M(3,3)*h*A];

InvAA = inv(AA);

F123 = zeros(6,1);
F = zeros(6,1);
  
for i = 2 : n

  F123(1:2) = A*Yy(:,i-1);
  F123(3:4) = A*Yy(:,i-1);
  F123(5:6) = A*Yy(:,i-1);
  
  F = InvAA * F123;
  
  Yy(:,i) = Yy(:,i-1) + h*((4/9 - sqrt(6)/36)*F(1:2) + (4/9 + sqrt(6)/36)*F(3:4) + 1/9*F(5:6));
  
end

Y(1,:) = Yy(1,:);

end

function Y = euler_(h, n)

Y = zeros(1, n);
y = [0; 1];
f = @(y) [y(2); -1/3 * y(2) - 1/4 * y(1)];

for i = 2 : n
  y = y + h * f(y);
  Y(i) = y(1);
end
end
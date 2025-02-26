% Hypertension/exercise model of vascular flow
% function hypertensionexercise
% SJP 21.viii.24
function hypertensionexercise

a1 = 1; a2b = 1; b1 = 0.01; b2 = 0.2; al = 1; pb = 20;

t = 0:0.01:60; % Time in years
pa = 1+(0.25*t/60);
h = 1+(0.001*t/60);
y0 = [0 0 1];
[t,y] = ode45(@(t,y) odefn(t,y,a1,a2b,b1,b2,al,pb,pa,h), t, y0);
s1 = y(:,1); s2 = y(:,2); k = y(:,3);

% Compute radius and flow
r = 0.1:0.1:2; % Range of r
RT = 0.2+(0.5./(r.^4))+0.3;
Pa0P0 = 1/(0.25+0.3); R = 0*pa;
f = (((r.^2)-1)/al)+log(pb+1);
for i = 1:length(t)
    fk = (exp(k(i)*f)-1)/pb;
    fr = (Pa0P0*(pa(i).*(0.3+(0.25./(r.^4)))./RT))-fk;
    R(i) = interp1(fr,r,0,'spline');
end
q = pa.*(R.^4);

set(0,'defaultAxesFontSize',16), set(0, 'DefaultLineLineWidth', 2)
figure(1)
subplot(2,1,1), plot(t, s1, 'r', t, s2, 'g'), ylabel('Response')
legend('s_1', 's_2', 'Location', 'northwest')
subplot(2,1,2), plot(t, pa, 'k', t, k, 'g', t, R, 'b', t, q, 'r')
xlabel('Time (years)'), ylabel('Response')
legend('ABP', 'Stiffness', 'Radius', 'Flow', 'Location', 'southwest')
end

% Differential function
function dydt = odefn(t,y,a1,a2b,b1,b2,al,pb,pa,h)
% Read variables
s1 = y(1); s2 = y(2); k = y(3);
pa = interp1(0:0.01:60,pa,t); h = interp1(0:0.01:60,h,t);
% Compute radius by interpolation, then pm
r = 0.1:0.1:2; % Range of r
RT = 0.2+(0.5./(r.^4))+0.3;
Pa0P0 = 1/(0.25+0.3);
f = (((r.^2)-1)/al)+log(pb+1);
fk = (exp(k*f)-1)/pb;
fr = (Pa0P0*(pa.*(0.3+(0.25./(r.^4)))./RT))-fk;
R = interp1(fr,r,0,'spline');
pm = (exp(k*((((R^2)-1)/al)+log(pb+1)))-1)/pb;
dydt(1) = a1*(pm-1);
dydt(2) = a2b*(h-1);
dydt(3) = -(b2*(k-1))+(b1*(s1-s2));
dydt = dydt(:);
end

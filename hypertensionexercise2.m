% Hypertension/exercise model of vascular flow
% function hypertensionexercise2
% SJP 23.viii.24
function hypertensionexercise2

a1 = 0.25; a2 = 0.25; b1 = 0.1; b2 = 1/3; b3 = 0.1; b4 = 1/3;
g1 = 1/3; g2 = 1/3; g3 = 1/3; al = 1; pb = 20;

t = 20:0.01:80; % Time in years
pa = 1+(0.25*(t-20)/60);
fE = 1-(0.5*(t-20)/60);
y0 = [0 0 1 1 1]; % s1 s2 F k e
[t,y] = ode45(@(t,y) odefn(t,y,a1,a2,b1,b2,b3,b4,g1,g2,g3,al,pb,pa,fE), t, y0);
s1 = y(:,1); s2 = y(:,2); F = y(:,3); kk = y(:,4); e = y(:,5); k = kk./e;

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
subplot(2,1,1), plot(t, s1, 'r', t, s2, 'g', t, F, 'k'), ylabel('Response')
legend('s_1', 's_2', 'Fitness', 'Location', 'northwest')
subplot(2,1,2), plot(t, pa, 'k', t, k, 'g', t, R, 'b', t, q, 'r', t, e, 'c')
xlabel('Age (years)'), ylabel('Response')
legend('ABP', 'Stiffness', 'Radius', 'Flow', 'Endothelium', 'Location', 'southwest')
end

% Differential function
function dydt = odefn(t,y,a1,a2,b1,b2,b3,b4,g1,g2,g3,al,pb,pa,fE)
% Read variables
s1 = y(1); s2 = y(2); F = y(3); kk = y(4); e = y(5); k = kk/e;
pa = interp1(20:0.01:80,pa,t); fE = interp1(20:0.01:80,fE,t);
% Compute radius by interpolation, then pm
r = 0.1:0.1:2; % Range of r
RT = 0.2+(0.5./(r.^4))+0.3;
Pa0P0 = 1/(0.25+0.3);
f = (((r.^2)-1)/al)+log(pb+1);
fk = (exp(k*f)-1)/pb;
fr = (Pa0P0*(pa.*(0.3+(0.25./(r.^4)))./RT))-fk;
R = interp1(fr,r,0,'spline');
pm = (exp(k*((((R^2)-1)/al)+log(pb+1)))-1)/pb;
dydt(1) = (a1*(pm-1))-(g1*s1);
dydt(2) = (a2*(F-1))-(g2*s2);
dydt(3) = -g3*(F-fE);
dydt(4) = -(b2*(kk-1))+(b1*(s1-s2));
dydt(5) = -(b4*(e-1))-(b3*(s1-s2));
dydt = dydt(:);
end

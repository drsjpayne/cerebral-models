% Hypertension/exercise model of vascular flow
% function [cbf] = hypertensionexercise3(dmax, ton, toff, he)
% SJP 28.viii.24
function [cbf] = hypertensionexercise3(dmax, ton, toff, he)

a0 = 1/60; a1 = 0.12; a2 = 0.12; b1 = 0.3; b2 = b1*a0/0.012; b3 = 0.25; b4 = 0.25; Fmin = 0.5;
g1 = 0.1; g2 = 0.1; g3 = 0.1; m12 = 0.28/0.2; al = 0.5; pb = 1.25; r1 = 0.2; r2 = 0.5; r3 = 1-r1-r2;

t = 20:0.1:80; % Time in years
if he == 1 % Hypertension only
    pa = 1+(dmax*0.5*(tanh((t-ton)/5)-tanh((t-toff)/5)));
    d = 0*dmax*0.5*(tanh((t-ton)/5)-tanh((t-toff)/5)); % Baseline value is 0.2
elseif he == 2 % Exercise only
    pa = 1+(0*dmax*0.5*(tanh((t-ton)/5)-tanh((t-toff)/5)));
    d = dmax*0.5*(tanh((t-ton)/5)-tanh((t-toff)/5)); % Baseline value is 0.2
elseif he == 3 % Combined
    pa = 1+(dmax(1)*0.5*(tanh((t-ton(1))/5)-tanh((t-toff(1))/5)));
    d = dmax(2)*0.5*(tanh((t-ton(2))/5)-tanh((t-toff(2))/5)); % Baseline value is 0.2
end
y0 = [0 0 0 1 0 1 1]; % state variables are s0 s1 s2 F f k1 e
[t,y] = ode15s(@(t,y) odefn(t,y,a0,a1,a2,b1,b2,b3,b4,g1,g2,g3,m12,Fmin,al,pb,r1,r2,r3,pa,d), t, y0);
s0 = y(:,1); s1 = y(:,2); s2 = y(:,3); F = y(:,4); k1 = y(:,6); e = y(:,7);
k = k1-(m12*e.*F.*d(:)); % arterial stiffness k = k1+k2

% Compute radius and flow
r = 0.1:0.1:2; % Range of r
RT = r1+(r2./(r.^4))+r3;
Pa0P0 = 1/((r2/2)+r3); R = 0*pa;
f = (((r.^2)-1)/al)+log(pb+1);
for i = 1:length(t)
    fk = (exp(k(i)*f)-1)/pb;
    fr = (Pa0P0*(pa(i).*(r3+((r2/2)./(r.^4)))./RT))-fk;
    R(i) = interp1(fr,r,0,'spline');
end
q = pa.*(R.^4); cbf = q;

set(0,'defaultAxesFontSize',16), set(0, 'DefaultLineLineWidth', 2)
figure(1)
subplot(2,2,1), plot(t, s0, 'y', t, s1, 'r', t, s2, 'g'), ylabel('Signals')
legend('s_0', 's_1', 's_2', 'Location', 'northwest'), title('(a)')
subplot(2,2,2), plot(t, k, 'k', t, e, 'b'), ylabel('Vascular health')
legend('Stiffness', 'Endothelium', 'Location', 'northwest'), title('(b)')
subplot(2,2,3), plot(t, F, 'k',t, d, 'b')
xlabel('Age (years)'), ylabel('Fitness and exercise'), axis([20 80 -0.01 1.01])
legend('Fitness', 'Exercise'), title('(c)')
subplot(2,2,4), plot(t, pa, 'k', t, R, 'b', t, q, 'r')
xlabel('Age (years)'), ylabel('Vascular parameters')
legend('ABP', 'Radius', 'Flow'), title('(d)')
end

% Differential function
function dydt = odefn(t,y,a0,a1,a2,b1,b2,b3,b4,g1,g2,g3,m12,Fmin,al,pb,r1,r2,r3,pa,d)
% Read variables
s0 = y(1); s1 = y(2); s2 = y(3); F = y(4); f = y(5); k1 = y(6); e = y(7);
pa = interp1(20:0.1:80,pa,t); d = interp1(20:0.1:80,d,t);
k = k1-(m12*e*F*d);
% Compute radius by interpolation, then pm etc.
r = 0.1:0.1:2; % Range of r
RT = r1+(r2./(r.^4))+r3;
Pa0P0 = 1/((r2/2)+r3);
FF = (((r.^2)-1)/al)+log(pb+1);
Fk = (exp(k*FF)-1)/pb;
Fr = (Pa0P0*(pa.*(r3+((r2/2)./(r.^4)))./RT))-Fk;
R = interp1(Fr,r,0,'spline');
pm = (exp(k*((((R^2)-1)/al)+log(pb+1)))-1)/pb;
dydt(1) = a0;
dydt(2) = max([(a1*(pm-1)) 0])-(g1*F*s1);
dydt(3) = (a2*e*d)-(g2*(1-F)*s2);
dydt(4) = -g3*(F-Fmin-f); if F > 1, dydt(4) = min([0 dydt(4)]); end
dydt(5) = d;
dydt(6) = (b1*(s0+s1-s2))-(b2*F*(k1-1));
dydt(7) = -(b3*(s1-s2))+(b4*F*(1-e)); if e > 1, dydt(7) = min([0 dydt(7)]); end
dydt = dydt(:);
end

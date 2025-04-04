clc
clear

%Constants

R = 8.314; 
g = 9.81;

%Initial conditions

T = 42e-3 % Thrust (N)
Pc = 150000; % chamber pressure (Pa)
Tc = 300; % chamber temperature (K)
MM = 0.028 ; % kg/mol
gamma = 1.4;
eps = 4; % area ratio

% Exit Pressure

expansion = @(Pe) ((gamma + 1) / 2)^(1 / (gamma - 1)) * ...
          (Pe / Pc)^(1 / gamma) * ...
          sqrt((gamma + 1) / (gamma - 1) * (1 - (Pe / Pc)^((gamma - 1) / gamma))) - (1/eps);

Pe = fsolve(expansion, 500) 

% alpha = 30;
% lambda = (1 + cos(alpha))/2;
% lambda=1;
% k = 1.4; 
% 
% 
% aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
% bb = 2/k;
% cc = (k-1)/k;
% 
% fun = @(x) x^bb*(1-x^cc)-aa;
% dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);
% 
% c1 = 0;
% c2 = 0.5;
% tol = 1e-6;
% err = tol + 1;
% it = 0;
% 
% while (it < 3 && err > tol ) 
%     it=it+1;
%     x = (c2+c1)/2; %stima dello zero
%     fc = fun(x);     
%     err=abs(fc); 
%     % scelta del nuovo estremo per l'eventuale ciclo successivo       
%     if (fc*fun(c1) > 0)
%           c1=x; 
%     else 
%           c2=x; 
%     end
% end
% xv = x;
% while (it< 5 && err> tol)
%    dfx = dfun(xv);
%    if dfx == 0
%       error(' Arresto per azzeramento di dfun');
%    else
%       xn = xv - fun(xv)/dfx;
%       err = abs(fun(xn));
%       it = it+1;
%       xv = xn;
%    end
% end
% 
% % Calcolo della pressione di uscita
% Pe = xn * Pc
% 
% Area Throat
At = (T/Pc)/(sqrt(2*(gamma^2/(gamma-1))*(2/(gamma+1))^((gamma+1)/(gamma-1)))*sqrt(1-(Pe/Pc)^((gamma-1)/gamma))+eps*(Pe-0)/Pc)


% Exit Velocity

v_exit = sqrt(2*gamma/(gamma-1)*R/MM*Tc*(1-(Pe/Pc)^((gamma-1)/gamma)))




% Exit Area

Ae = eps*At
re=sqrt(Ae/pi)
De=2*re
% Mass Flow Rate

m_dot=(T-Ae*(Pe-0))/v_exit


% Isp

Isp = T/m_dot/g


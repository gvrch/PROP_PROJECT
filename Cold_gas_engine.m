clc
clear

%Constants

R = 8.314; 
g = 9.81;

%Initial conditions

T = 42e-3; % Thrust (N)
Pc = 100000; % chamber pressure (Pa)
Tc = 300; % chamber temperature (K)
MM = 0.028 ; % kg/mol
gamma = 1.4;
eps = 4; % area ratio

% Exit Pressure


expansion = @(Pe) ((gamma + 1) / 2)^(1 / (gamma - 1)) * ...
          (Pe / Pc)^(1 / gamma) * ...
          sqrt((gamma + 1) / (gamma - 1) * (1 - (Pe / Pc)^((gamma - 1) / gamma))) - (1/eps);

options = optimset('TolX',1e-10);
Pe = fsolve(expansion, 10000, options) 

% Area Throat
At = (T/Pc)/(sqrt(2*(gamma^2/(gamma-1))*(2/(gamma+1))^((gamma+1)/(gamma-1)))*sqrt(1-(Pe/Pc)^((gamma-1)/gamma))+eps*(Pe-0)/Pc)
rt=sqrt(At/pi)
Dt=2*rt

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

c_star = At*Pc/m_dot

% to-do:
% nitrogen storage information, data on thruster everything
% uncertainties: nozzle features (diameter, exit area), Pc
% losses: nozzle, pressure
% monte carlo
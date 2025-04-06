%% Vaporizing Liquid Microthrusters V4
%
% -------------------------------------------------------------------------
clear
clc
close all

%% DATA
eps     = 9;            % Rapporto tra le aree (ugello convergente-divergente)  [adim]
k       = 1.327;        % Rapporto dei calori specifici                         [adim]
g       = 9.81;         % Gravity acceleration                                  [m/s^2]
Pc      = 230000;       % Evaporation chamber Pressure                          [Pa]
R       = 8.314;        % Gas Constant                                          [J/mol*K]
Mmol    = 0.018;        % Water Molar Mass                                      [kg/mol]
rho_l   = 1000;         % Density liquid water                                  [kg/m^3]
rho_v   = 0.6;          % Density gas water                                     [kg/m^3]
cp_l    = 4.1838e3;     % Specific Heat liquid water                            [J/kgK]
cp_v    = 2.0256e3;     % Specific Heat gas water                               [J/kgK]
lambda  = 2260e+3;      % Water latent evaporation heat                         [J/kg]
T_eb    = 124 + 273;    % Temperature of ebollition                             [K]
T_inj   = 25 + 273.15;  % temperature of Injection                              [K]
alpha   = deg2rad(30);  %                                                       [deg]
beta    = deg2rad(45);  %                                                       [deg]
Cd      = 0.7;          % Drag coefficient                                      [adim]
dP_inj  = 0.1;          % Pressure lost for injection   !!!CHECK!!!
T       = 5e-3;         % Wanted Thrust                                         [N]
dv      = 0.5;          % Wanted Delta v                                        [m/s]
mass    = 4;            % S/C mass                                              [kg]
N       = 100;          % Number of sample for Q_dot                            [adim]


%% Calculate the Pressure Ratio (Pe/Pc)

% equation
aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
bb = 2/k;
cc = (k-1)/k;

% Resolve with Newton method
fun  = @(x) x^bb*(1-x^cc)-aa;
dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);

c1  = 0;
c2  = 0.5;
tol = 1e-6;
err = tol + 1;
it  = 0;

while (it < 3 && err > tol ) 
    it  = it+1;
    x   = (c2+c1)/2; %stima dello zero
    fc  = fun(x);     
    err =abs(fc); 
    % scelta del nuovo estremo per l'eventuale ciclo successivo       
    if (fc*fun(c1) > 0)
          c1=x; 
    else 
          c2=x; 
    end
end
xv = x;
while (it< 10 && err> tol)
   dfx = dfun(xv);
   if dfx == 0
      error(' Arresto per azzeramento di dfun');
   else
      xn = xv - fun(xv)/dfx;
      err = abs(fun(xn));
      it = it+1;
      xv = xn;
   end
end

% Compute exit pressure value
Pe = xn * Pc;

%% Compute throat and exit area 

CT = k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k)))+Pe/Pc*eps;
At = T/Pc/CT;
r_t = sqrt(At/pi);
Ae = eps*At;
r_e = sqrt(Ae/pi);

%% Dimensionate nozzle
eps_conv = 1/0.2 * sqrt((2/(k+1)*(1+(k-1)/2*0.2^2))^((k+1)/(k-1)));
Ac = eps_conv * At;
r_c = sqrt(Ac/pi);

L_conv = 1/2*(2*(r_c-r_t)/tan(beta)); % ipotesi che Ae è uagule a quwlla di ingresso nel convergente
L_div = 1/2*(2*(r_e-r_t)/tan(alpha));
L_nozzle = L_conv + L_div;
l = 1/2*(1+cos(alpha));

%% Study the flow in function of Q_dot

% Define some usefull constant value 
K       = sqrt(2*k/(k-1) * R/Mmol * (1 - (Pe/Pc)^((k-1)/k))); % V_e = k * sqrt(T_c)
aa      = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v;
cc      = cp_v*((T-Pe*Ae)/K)^2;
Q_min   = sqrt(4*aa*cc);

m_dot   = zeros(N,2);
T_c     = zeros(N,2);
v_e     = zeros(N,2);
T_e     = zeros(N,2);
Isp     = zeros(N,2);
A_inj   = zeros(N,2);
Q_dot   = linspace(Q_min,13,N);

for i = 1:length(Q_dot)
    
    % mass flow rate
    m_dot(i,1) = (Q_dot(i) + sqrt(Q_dot(i)^2 - Q_min^2))/(2*aa);
    m_dot(i,2) = (Q_dot(i) - sqrt(Q_dot(i)^2 - Q_min^2))/(2*aa);
    
    % Combustion Chamber Temperature
    T_c(i,1) = ((T-Pe*Ae)/(m_dot(i,1) *K))^2;
    T_c(i,2) = ((T-Pe*Ae)/(m_dot(i,2) *K))^2;
    
    % Exit Velocity
    v_e(i,1) = K *sqrt(T_c(i,1));
    v_e(i,2) = K *sqrt(T_c(i,2));

    % Exit Temperature
    T_e(i,1) = T_c(i,1)*(Pe/Pc)^((k-1)/k);
    T_e(i,2) = T_c(i,2)*(Pe/Pc)^((k-1)/k);
    
    % Specific Impulse
    Isp(i,1) = T/(m_dot(i,1)*g);
    Isp(i,2) = T/(m_dot(i,1)*g);

    % Injection Area
    A_inj(i,1) = m_dot(i,1)/(Cd*2*sqrt(dP_inj*Pc*rho_l));
    A_inj(i,2) = m_dot(i,2)/(Cd*2*sqrt(dP_inj*Pc*rho_l));
end

%% Compute specific impulse


%% PLOTS

figure

% Mass flow rate 
subplot(2,4,1)
plot(Q_dot,m_dot(:,1), "LineWidth",1)
hold on
plot(Q_dot,m_dot(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Mass flow Rate [kg/s]');
legend('Solution #1','Solution #2');

% Combustion Chamber Temperature
subplot(2,4,2)
plot(Q_dot,T_c(:,1), "LineWidth",1)
hold on
plot(Q_dot,T_c(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Combustion Chamber Temperature [K]');
legend('Solution #1','Solution #2');

% Exit Temperature
subplot(2,4,3)
plot(Q_dot,T_e(:,1), "LineWidth",1)
hold on
plot(Q_dot,T_e(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Exit Temperature [K]');
legend('Solution #1','Solution #2');

% Exit Velocity
subplot(2,4,4)
plot(Q_dot,v_e(:,1), "LineWidth",1)
hold on
plot(Q_dot,v_e(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Exit Velocity [m/s]');
legend('Solution #1','Solution #2');

% Specific Impulse
subplot(2,4,5)
plot(Q_dot,Isp(:,1), "LineWidth",1)
hold on
plot(Q_dot,Isp(:,2),"LineStyle","--", "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Specific Impulse [s]');
legend('Solution #1','Solution #2');

% Injection Area
subplot(2,4,6)
plot(Q_dot,A_inj(:,1), "LineWidth",1)
hold on
plot(Q_dot,A_inj(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Injection Area [m^2]');
legend('Solution #1','Solution #2');

% Nozzle
figure


%% VAPORIZING LIQUID MICROTHRUSTER
% Margi e Ludo 

% CHECK: possibili migliorie
% - k del vapore acqueo per valore più adatto a Tc(Q_dot)
% - densità del vapore d'acqua rho_v(T)

%%

clc
clear
close all

%% NOMINAL OPERATION

% Dati del problema
eps = 10;        % Rapporto tra le aree (ugello convergente-divergente)
k = 1.327;       % Rapporto dei calori specifici
Pc = 230000;     % Pressione in camera di combustione (Pa)
R = 8.314;       % Costante specifica del gas (J/kg*K), da definire correttamente per il propellente
Mmol = 0.018;    % Massa molare del gas (kg/mol), da definire correttamente
rho_l = 1000;             % kg/m^3
rho_v = 0.6;                
% Q_dot = 20;                % W
T = 5e-3;                   % N
cp_l = 4186;                % J/kgK
cp_v = 1910;                % J/kgK
lambda = 2260e+3;           % J/kg
T_eb = 124 + 273.15;        % K
T_inj = 25 + 273.15;        % K 
alpha = deg2rad(30);
beta = deg2rad(45);
Cd = 0.7;                   % 
dv = 0.5;                   % m/s
mass = 4;                   % kg

% Funzione anonima per risolvere il rapporto tra le pressioni
f = @(x) -1/eps + ((k+1)/2)^(1/(k-1)) * x^(1/k) * sqrt((k+1)/(k-1) * (1-x^((k-1)/k)));

% Risoluzione per il rapporto di pressione x = Pe/Pc
x = fsolve(f, 0.001);  % Stima iniziale bassa (tipico Pe/Pc è molto piccolo)

% Calcolo della pressione di uscita
Pe = x * Pc;

% Calcolo di area di gola e area di uscita
CT = k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k)))+Pe/Pc*eps;
At = T/Pc/CT;
r_t = sqrt(At/pi);
Ae = eps*At;
r_e = sqrt(Ae/pi);

eps_conv = 1/0.2 * sqrt((2/(k+1)*(1+(k-1)/2*0.2^2))^((k+1)/(k-1)));
Ac = eps_conv * At;
r_c = sqrt(Ac/pi);

L_conv = 1/2*(2*(r_c-r_t)/tan(beta));
L_div = 1/2*(2*(r_e-r_t)/tan(alpha));
L_nozzle = L_conv + L_div;
l = 1/2*(1+cos(alpha));

Tc_v = [T_eb+10:10:580.15+273.15]; % K   % 273.15 è la T operazionale max del Silicio

TE = [];
Q = [];
M = [];
V = [];
ISP = [];
A_INJ = [];
Mq = [];

for Tc = Tc_v
    % Q_dot = 10;
    % m_dot_q = Q_dot/(cp_l*(T_eb - T_inj) + lambda + cp_v*(Tc - T_eb));
    ve = sqrt(2*k/(k-1)*R/Mmol*Tc*(1-(Pe/Pc)^((k-1)/k)));
    m_dot = (T-Pe*Ae)/(ve*l);
    Q_dot = m_dot*(cp_l*(T_eb - T_inj) + lambda + cp_v*(Tc - T_eb));
    Isp = T/(m_dot*9.81);
    A_inj = m_dot/(Cd*2*sqrt(0.10*Pc*rho_l));
    
    dt = dv * mass / T;
    m_fuel_1 = m_dot*dt;
    v_fuel = m_fuel_1/rho_l;
    r_tank_sph = (3/4/pi*v_fuel)^(1/3);
    
    Te = Tc * (Pe/Pc)^((k-1)/k);
    
    Q = [Q, Q_dot];
    M = [M, m_dot];
    %Mq = [Mq, m_dot_q];
    V = [V, ve]; 
    ISP = [ISP, Isp];
    A_INJ = [A_INJ, A_inj];
    TE = [TE, Te];
end


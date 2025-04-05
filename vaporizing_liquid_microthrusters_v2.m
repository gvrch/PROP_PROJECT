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
Pc = 200000;     % Pressione in camera di combustione (Pa)
R = 8.314;       % Costante specifica del gas (J/kg*K), da definire correttamente per il propellente
Mmol = 0.028;    % Massa molare del gas (kg/mol), da definire correttamente
rho_l = 1000;             % kg/m^3
rho_v = 0.6;                
%Q_dot = 20;                % W
T = 5e-3;                   % N
cp_l = 4186;                % J/kgK
cp_v = 1910;                % J/kgK
lambda = 2260e+3;           % J/kg
T_eb = 124 + 273;        % K
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

L_conv = 1/2*(2*(r_c-r_t)/tan(beta)); % ipotesi che Ae è uagule a quwlla di ingresso nel convergente
L_div = 1/2*(2*(r_e-r_t)/tan(alpha));
L_nozzle = L_conv + L_div;
l = 1/2*(1+cos(alpha));

TC = [];
M = [];
V = [];
ISP = [];
A_INJ = [];
TE = [];

for Q_dot = [1:1:20]
% Calcolo della velocità di scarico
K = sqrt(2*k/(k-1)*R/Mmol*(1-(Pe/Pc)^((k-1)/k)));
m_dot = Q_dot / (cp_l*(T_eb - T_inj) + lambda + cp_v*((T-Pe*Ae)/(l*K))^2 - cp_v*T_eb);
Tc = ((T-Pe*Ae)/(m_dot*K*l))^2;
ve = K*sqrt(Tc);

Isp = T/(m_dot*9.81);
A_inj = m_dot/(Cd*2*sqrt(0.10*Pc*rho_l));

dt = dv * mass / T;
m_fuel_1 = m_dot*dt;
v_fuel = m_fuel_1/rho_l;
% m_final = mass/exp(dv/(Isp*9.81));
% m_fuel_2 = mass - m_final;

Te = Tc * (Pe/Pc)^((k-1)/k)

TC = [TC, Tc];
M = [M, m_dot];
V = [V, ve]; 
ISP = [ISP, Isp];
A_INJ = [A_INJ, A_inj];
TE = [TE, Te];
end
TC = TC - 273.15;
TE = TE - 273.15;
% cp_l = 4186;                %J/kgK
% cp_v = 1910;                %J/kgK
% lambda = 2260e+3;           %J/kg
% 
% H_v = 3488.42e+3; % A 500 C
% H_l = 83.378e+3; 
% 
% m_dot = 5e-6; %kg/s
% Q_dot = 17; %W
% T_imposta = 500 + 273.15; %W
% 
% T_eb = 100.3 + 273; % K
% T_inj = 25 + 273.15; % K 
% 
% 
% Tc_fun = (Q_dot-m_dot*cp_l*(T_eb-T_inj)-m_dot*lambda)/(m_dot*cp_v) + T_eb
% Q_fun = m_dot*cp_l*(T_eb - T_inj) + m_dot*lambda + m_dot*cp_v*(T_imposta - T_eb)
% 
% Q_nist = m_dot*(H_v-H_l)

%% MONTECARLO

%% Definizione parametri
k_p     = 1;        % Coefficiente di perdita tubi
k_v     = 1.3;      % Coefficiente di perdita valvola
k_inj   = 1.17;     % Coefficiente di perdita iniettore
% T_c    = 300;     % Chamber Temp in K
g       = 9.81;     % Accelerazione di gravità
rho     = rho_l;    % Densità del propellente liquido (kg/m^3)

P_cc_nom    = Pc;           % Pressione nominale in camera di combustione (Pa)
A_t_nom     = At;           % Area ugello (m^2)
A_inj_nom   = A_inj;        % Area nominale iniettore (m^2)
m_dot_nom   = M(10);        % Flusso massico nominale (kg/s)
Q_dot_nom   = 10;           %Flusso di calore nominale (W)
A_tube      = 0.00042;      % Area del tubo di alimentazione (m^2) CAMBIARE
A_v         = 32e-5;        % Area delle valvole (m^2)
tol         = 1e-2;         % Tolleranza convergenza


c_star = P_cc_nom * A_t_nom / m_dot_nom;   % Velocità caratteristica del propellente (m/s)
%sqrt(R*T_c/k*((k+1)/2)^((k+1)/(k-1)));

%% Definizione nuovo punto di lavoro

P_cc = P_cc_nom;
P_tank = P_cc_nom*1.2;
m_dot2 = 15;
m_dot1 = 28;

N = 10000;
rt_avg = sqrt(A_t_nom/pi);
rt_std = 1e-6;
rt_val = normrnd(rt_avg, rt_std, [N,1]);
At_val = rt_val.^2.*pi;

r_inj_std = 1e-6;       
r_inj_avg = sqrt(A_inj_nom/pi); 
r_inj_val = normrnd(r_inj_avg, r_inj_std, [N,1]);
A_inj_val = r_inj_val.^2.*pi;

Q_dot_min = Q_dot_nom*0.95;       
Q_dot_max = Q_dot_nom*1.01;
Q_dot_val = Q_min + (Q_max - Q_min) * rand(1, N);


K_tot_val = k_v*A_v/A_tube + k_p + k_inj.*A_inj_val/A_tube;

T_MC = [];
for ii = 1:N
    
    A_t  = At_val(ii);
    k_in = K_tot_val(ii);
    Q_dot = Q_dot_val(ii);
    it = 0;

    P_cc = P_cc_nom;
    m_dot2 = 5;
    m_dot1 = 19;

    while(abs(m_dot2 - m_dot1) > tol && it < 1000)
        it = it + 1;

        m_dot1 =  P_cc*A_t/c_star;
        v  = m_dot1 /( A_tube* rho);
        deltaP= 0.5*rho*v^2*k_in;
        P_cc = P_tank - deltaP;
        
        m_dot2= P_cc*A_t/c_star;
       
        if it > 999
            disp('It did not converge')
        end   

    end

    T_c = (Q_dot - m_dot1*cp_l*(T_eb - T_inj) - m_dot*lambda)/(m_dot*cp_v) + T_eb; 

    % calcolo Pe
    aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
    bb = 2/k;
    cc = (k-1)/k;

    fun = @(x) x^bb*(1-x^cc)-aa;
    dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);

    c1 = 0;
    c2 = 0.5;
    toll = 1e-6;
    err = toll + 1;
    it = 0;
    iter = 0;
    while (iter < 3 && err > toll ) 
        iter=iter+1;
        x = (c2+c1)/2; %stima dello zero
        fc = fun(x);     
        err=abs(fc); 
        % scelta del nuovo estremo per l'eventuale ciclo successivo       
        if (fc*fun(c1) > 0)
              c1=x; 
        else 
              c2=x; 
        end
    end
    xv = x;
    iter = 0;
    while (iter< 5 && err> toll)
       dfx = dfun(xv);
       if dfx == 0
          error(' Arresto per azzeramento di dfun');
       else
          xn = xv - fun(xv)/dfx;
          err = abs(fun(xn));
          iter = iter+1;
          xv = xn;
       end
    end
    Pe = xv*P_cc;
    %
    v_e = sqrt(2*k/(k-1)*R/Mmol*T_c*(1-(Pe/P_cc)^((k-1)/k)));
    T_MC(ii) =  m_dot1*v_e + Ae*Pe;
    Isp_MC(ii) = T_MC(ii)/(m_dot2*g);
    if ii > 1
        T_avg(ii) = T_avg(ii-1) + 1/ii*(T_MC(ii) - T_avg(ii-1));
        Isp_avg(ii) = Isp_avg(ii-1) + 1/ii*(Isp_MC(ii) - Isp_avg(ii-1));
    else
        T_avg(ii) = T_MC(ii);
        Isp_avg(ii) = Isp_MC(ii);
    end
    T_std(ii) = std(T(1:ii));
    Isp_std(ii) = std(Isp(1:ii));
end

subplot(2,2,1)
plot(T_avg)
title T_{avg}
subplot(2,2,2)
plot(T_std)
title T_{std}
subplot(2,2,3)
plot(Isp_avg)
title Isp_{avg}
subplot(2,2,4)
plot(Isp_std)
title Isp_{std}

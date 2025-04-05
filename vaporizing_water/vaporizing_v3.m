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
eps = 100;              % Rapporto tra le aree (ugello convergente-divergente)
k = 1.327;              % Rapporto dei calori specifici
Pc = 200000;            % Pressione in camera di combustione (Pa)
R = 8.314;              % Costante specifica del gas (J/mol*K), da definire correttamente per il propellente
Mmol = 0.018;           % Massa molare del gas (kg/mol), da definire correttamente
rho_l = 1000;           % kg/m^3
rho_v = 0.6;                
%Q_dot = 20;            % W
T = 5e-3;               % N
cp_l = 4.1838e3;        % J/kgK
cp_v = 2.0256e3;        % J/kgK
lambda = 2260e+3;       % J/kg
T_eb = 124 + 273;       % K
T_inj = 25 + 273.15;    % K 
alpha = deg2rad(30);
beta = deg2rad(45);
Cd = 0.7;               % 
dv = 0.5;               % m/s
mass = 4;               % kg

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
 l= 1;
% l = 1/2*(1+cos(alpha));

syms v_e_var m_dot_var Pc_var R_var MM_var k_var T_var positive real
syms Pe_var Ae_var Q_dot_var C_pl_var C_pv_var T_eb_var T_inj_var T_c_var lambda_var positive real

eq(1) = T_var == m_dot_var * v_e_var + Pe_var * Ae_var;
eq(2) = Q_dot_var == m_dot_var * (C_pl_var * (T_eb_var - T_inj_var) + lambda_var + C_pv_var*(T_c_var - T_eb_var));
eq(3) = v_e_var == sqrt((2*k_var/(k_var-1))*R_var/MM_var*T_c_var*(1-(Pe_var/Pc_var)^((k_var-1)/k_var)));

eq = subs(eq,[Pc_var,R_var,MM_var,k_var, T_var,Pe_var,Ae_var,C_pl_var,C_pv_var,T_eb_var,T_inj_var,lambda_var],[Pc,R,Mmol,k,T,Pe,Ae,cp_l,cp_v,T_eb,T_inj,lambda]);


TC = [];
M = [];
VE = [];
ISP = [];
A_INJ = [];
TE = [];


for Q_dot = 1:1:20

    eq_temp = subs(eq,Q_dot_var,Q_dot);
    %sol = solve(eq_temp,[T_c_var, v_e_var, m_dot_var]);
         sol = vpasolve(eq_temp,[T_c_var, v_e_var, m_dot_var], [100 800; 200 1500; 1e-8 1e-4]);
    %     [250 800; 400 1500; 1e-8 1e-4]
if isempty(sol.T_c_var)
    TC=[TC,0];
end
if isempty(sol.m_dot_var)
    M=[M,0];
end
if isempty(sol.v_e_var)
    VE=[VE,0];
end
TC = [TC,double(sol.T_c_var)];
M  = [M,double(sol.m_dot_var)];
VE = [VE,double(sol.v_e_var)];  

end

TC
M
VE


% 
% K1 = sqrt(2*k/(k-1)* R / Mmol * (1-(Pe/Pc)^((k-1)/k))); % First constant (K)
% K2 = cp_l * (T_eb - T_inj) + lambda - cp_v*T_eb;        % Second constant (N)
% K3 = cp_v * ((T-Pe*Ae)/(l*K1))^2;                        % Third constant (Z)
%  % m1 = Q_dot+sqrt(Q_dot^2-4*K2*K3)/(2*K2);
%     % m2 = Q_dot-sqrt(Q_dot^2-4*K2*K3)/(2*K2);
%     fun =@(m_dot) - Q_dot * m_dot + m_dot^2 * K2 + K3;
%     %fplot(fun,[0,1e-4]);
%     dfun = @(m_dot) - Q_dot + 2*m_dot * K2 ;
% 
%     c1 = 0;
%     c2 = 0.5;
%     toll = 1e-8;
%     err = toll + 1;
%     it = 0;
%     iter = 0;
%     while (iter < 3 && err > toll ) 
%         iter=iter+1;
%         x = (c2+c1)/2; %stima dello zero
%         fc = fun(x);     
%         err=abs(fc); 
%         % scelta del nuovo estremo per l'eventuale ciclo successivo       
%         if (fc*fun(c1) > 0)
%               c1=x; 
%         else 
%               c2=x; 
%         end
%     end
%     xv = x;
%     iter = 0;
%     while (iter< 8 && err> toll)
%        dfx = dfun(xv);
%        if dfx == 0
%           error(' Arresto per azzeramento di dfun');
%        else
%           xn = xv - fun(xv)/dfx;
%           err = abs(fun(xn));
%           iter = iter+1;
%           xv = xn;
%        end
%    end
% 
% 
% M1 = [M1, xv];
% %M2 = [M2, m2];

% m_dot = fsolve (f1,1e-6, options);

% syms m_dot_var real
% eq = Q_dot == m_dot_var * K2 + 1/m_dot_var * K3;
% m_dot_sol = double(solve(eq,m_dot_var));

%%
clear all
close all
clc

%%
addpath 'CEAM 1.0.0.0'

%%
% x=CEA('problem','hp','equilibrium','o/f', 6,'case','rocket','p,bar', 200*1e-5,...
% 'reactants',...
% 'fuel','H2(L)', 'H', 2, 'h,kj/kg', 5000,'wt%',100, 't(k)', 300, ...
% 'oxid','O2(L)', 'O', 2, 'h,kj/kg', 4000,'wt%',100,'t(k)',300,'output','massf', 'transport','end');

%% Definizione parametri
k_p = 1;    % Coefficiente di perdita primaria
k_v = 1.3;    % Coefficiente di perdita valvola
k_in_nom = 0.6;   % Coefficiente di ingresso
T_c = 300;      % Chamber Temp in K
gamma = 1.4; 
R = 8.314; 
g = 9.81;
MM = 0.028;



Ae = 1.6999e-06;
Pe = 2.9495e+03;


% m_dot_nom = 30; % Portata nominale (kg/s)

P_cc_nom = 150000; % Pressione nominale in camera di combustione (Pa)
A_t_nom = 1.6999e-07; % Area ugello (m^2)
A_tube = 0.00042; % Area del tubo di alimentazione (m^2)
rho = 807; % Densità del propellente liquido (kg/m^3)
tol = 1e-2; % Tolleranza convergenza

%%
m_dot_nom   = 6.1722e-05;
v_prel      = m_dot_nom /(A_tube*rho);
deltaP_prel = 0.5*rho*v_prel^2*(k_in_nom + k_p + k_v);

c_star = 435.8794; %sqrt(R*T_c/gamma*((gamma+1)/2)^((gamma+1)/(gamma-1))); % Velocità caratteristica del propellente (m/s)


%% Definizione nuovo punto di lavoro


P_cc = P_cc_nom;
P_tank = P_cc_nom + deltaP_prel;
m_dot2 = 5;
m_dot1 = 19;

A_t_vect = 0.99*A_t_nom:0.001*A_t_nom:1.01*A_t_nom;
K_in_vect = 0.99*k_in_nom:0.001*k_in_nom:1.01*k_in_nom;
T = [];
for A_t = A_t_vect

    for k_in = K_in_vect
    it = 0;

    P_cc = P_cc_nom;
    m_dot2 = 5;
    m_dot1 = 19;

    while(abs(m_dot2 - m_dot1) > tol && it < 1000)
        it = it + 1;

        m_dot1 =  P_cc*A_t/c_star;
        v  = m_dot1 /( A_tube* rho);
        deltaP= 0.5*rho*v^2*(k_in + k_p + k_v);
        P_cc = P_tank - deltaP;
        
        m_dot2= P_cc*A_t/c_star ;
       
        if it > 999
            disp('It did not converge')
        end   

    end
    v_e = sqrt(2*gamma/(gamma-1)*R/MM*T_c*(1-(Pe/P_cc)^((gamma-1)/gamma)));
    T = [T , m_dot1*v_e + Ae*Pe];
   
    end
end

T_mean  = mean(T)
T_sigma = std(T) 
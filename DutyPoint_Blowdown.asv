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
k_in = 0.6;   % Coefficiente di ingresso
c_star = 1600; % Velocità caratteristica del propellente (m/s)
% m_dot_nom = 30; % Portata nominale (kg/s)

P_cc_nom = 14e5; % Pressione nominale in camera di combustione (Pa)
A_t = 0.016; % Area ugello (m^2)
A_tube_ox = 0.00042; % Area del tubo di alimentazione (m^2)
% A_tube_f = 0.00030; % Area del tubo di alimentazione (m^2)
rho_ox = 1000; % Densità del propellente liquido (kg/m^3)
rho_f = 1200;
tol = 1e-2; % Tolleranza convergenza

%% o/f
OF = 2.6;
A_tube_f = A_tube_ox*sqrt(rho_ox/rho_f)/OF;

%% Definizione nuovo punto di lavoro

P_tank = 30e+5; % Nuova pressione nel serbatoio (Pa)
P_cc = P_cc_nom;
m_dot2 = 5;
m_dot1 = 19;

deltaP_vec = [];
mdot_vec = [];
m_ox = [];
m_f = [];


for P_tank = 30e+5:-100:20e+5
    m_dot2 = 5;
    m_dot1 = 19;
    it = 0;
    while(abs(m_dot2 - m_dot1) > tol && it < 1000)
        it = it + 1;
        deltaP = P_tank - P_cc;

        v_ox = sqrt(2 * deltaP / rho_ox / (k_in + k_p + k_v));
        v_f = sqrt(2 * deltaP / rho_f / (k_in + k_p + k_v));
        m_dot1 = rho_ox * v_ox * A_tube_ox + rho_f * v_f * A_tube_f;
    
        P_cc = c_star * m_dot1 / A_t;
        deltaP = P_tank - P_cc;
        v_ox = sqrt(2 * deltaP / rho_ox / (k_in + k_p + k_v));
        v_f = sqrt(2 * deltaP / rho_f / (k_in + k_p + k_v));
        m_dot2 = rho_ox * v_ox * A_tube_ox + rho_f * v_f * A_tube_f;
        
        if it > 999
            disp('Il calcolo non converge')
        end   

    end
         % disp(v_f)
         % disp(v_ox)
         m_ox = [m_ox, rho_ox * v_ox * A_tube_ox];
         m_f =  [m_f, rho_f * v_f * A_tube_f];
         deltaP_vec = [deltaP_vec, deltaP];
         mdot_vec = [mdot_vec, m_dot2];
% disp(P_tank)
% disp(['Punto di lavoro trovato dopo ', num2str(it), ' iterazioni.'])
% disp(['Pressione camera di combustione: ', num2str(P_cc / 1e5), ' bar'])
% disp(['Portata massica: ', num2str(m_dot2), ' kg/s'])
end


figure
scatter(mdot_vec, deltaP_vec)

figure
plot(mdot_vec, deltaP_vec)
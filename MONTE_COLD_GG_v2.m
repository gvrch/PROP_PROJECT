%%
clear
close all
clc

%% Definizione parametri
k_p = 1;        % Coefficiente di perdita primaria
k_v = 1.3;      % Coefficiente di perdita valvola
k_inj = 0.6;    % Coefficiente di ingresso
T_c = 300;      % Chamber Temp in K
gamma = 1.4; 
R = 8.314; 
g = 9.81;
MM = 0.028;



Ae = 1.6999e-06;
Pe = 2.9495e+03;


P_cc_nom = 150000;          % Pressione nominale in camera di combustione (Pa)
A_t_nom = 1.6999e-07;       % Area ugello (m^2)
A_tube = 0.00042;           % Area del tubo di alimentazione (m^2)
A_v = 32e-5;  
rho = 807;                  % Densità del propellente liquido (kg/m^3)
tol = 1e-2;                 % Tolleranza convergenza

c_star = 435.8794;          %sqrt(R*T_c/gamma*((gamma+1)/2)^((gamma+1)/(gamma-1))); % Velocità caratteristica del propellente (m/s)


%% Definizione nuovo punto di lavoro


P_cc = P_cc_nom;
P_tank = P_cc_nom*1.2;
m_dot2 = 5;
m_dot1 = 19;

N = 10000;
rt_avg = sqrt(A_t_nom/pi);
rt_std = 1e-6;
rt_val = normrnd(rt_avg, rt_std, [N,1]);
At_val = rt_val.^2.*pi;


r_inj_std = 1e-6;       
r_inj_avg = 500e-6; 
r_inj_val = normrnd(r_inj_avg, r_inj_std, [N,1]);
A_inj_val = r_inj_val.^2.*pi;

K_tot_val = k_v*A_v/A_tube + k_p + k_inj.*A_inj_val/A_tube;


T = [];
for ii = 1:N
    
    A_t  = At_val(ii);
    k_in = K_tot_val(ii);
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
    T(ii) =  m_dot1*v_e + Ae*Pe;
    Isp(ii) = T(ii)/(m_dot2*g);
    if ii > 1
        T_avg(ii) = T_avg(ii-1) + 1/ii*(T(ii) - T_avg(ii-1));
        Isp_avg(ii) = Isp_avg(ii-1) + 1/ii*(Isp(ii) - Isp_avg(ii-1));
    else
        T_avg(ii) = T(ii);
        Isp_avg(ii) = Isp(ii);
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
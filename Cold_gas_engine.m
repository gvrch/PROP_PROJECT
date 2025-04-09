clc
clear
close all

%Constants

R = 8.314; 
g = 9.81;


%Initial conditions

T = 5e-3;           % Thrust (N)
Pc = 150000;        % chamber pressure (Pa) 1bar-4bar
Tc = 300;           % chamber temperature (K) fixed
MM = 0.028 ;        % kg/mol
k = 1.4;
eps = 10;            % area ratio from 4-20
DeltaV = 0.5;       % requested total DV (m/s)
Mass = 4;           % 3U cubesat mass
rho = 1.176;
t_burn = DeltaV/(T/Mass);


RN = R/MM;

% feeding system
k_p = 1.1;  % P drop coeff of the pipes
k_v = 1.3;  % P drop valves
k_inj = 0.6;    % P_drop injector

% from paper but to be fixed
A_p = 1e-3^2*pi;
A_v = 0.5e-3^2*pi;
r_inj = 250e-6;
A_inj = r_inj^2*pi;

%% Nominal Design


aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
bb = 2/k;
cc = (k-1)/k;

fun = @(x) x^bb*(1-x^cc)-aa;
dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);

c1 = 0;
c2 = 0.5;
tol = 1e-6;
err = tol + 1;
it = 0;

while (it < 3 && err > tol ) 
    it=it+1;
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
while (it< 5 && err> tol)
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

% Calcolo della pressione di uscita
Pe = xn * Pc;

% temp in exit
Te = Tc*(Pe/Pc)^((k-1)/k);

% Area Throat
At = (T/Pc)/(sqrt(2*(k^2/(k-1))*(2/(k+1))^((k+1)/(k-1)))*sqrt(1-(Pe/Pc)^((k-1)/k))+eps*(Pe-0)/Pc);
rt=sqrt(At/pi);
Dt=2*rt;

% Exit Velocity
v_exit = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pe/Pc)^((k-1)/k)));
M_exit = v_exit/sqrt(k*RN*Te)


% Exit Area
Ae = eps*At;
re=sqrt(Ae/pi);
De=2*re;

% Mass Flow Rate
m_dot_nom=(T-Ae*(Pe-0))/v_exit;

% Isp
Isp_nom = T/m_dot_nom/g;

%cstar
c_star = At*Pc/m_dot_nom;

%% Preliminary sizing
% mass
MR_nom = exp(DeltaV/(Isp_nom*g));
M_N_used_nom = Mass*(MR_nom-1);

% pressure losses
v_p_nom = m_dot_nom/(rho*A_p);
Delta_P_tot_nom = 0.5*rho*v_p_nom^2*(k_p + A_p/A_v*k_v + A_p/A_inj*k_inj);
Delta_P_plen_nom = 0.5*rho*v_p_nom^2*A_p/A_inj*k_inj;
P_plen_nom = Pc + Delta_P_plen_nom;
P_f_nom = Pc + Delta_P_tot_nom;
P_f_nom = P_f_nom*1.2;

V_plen_nom = M_N_used_nom*RN*Tc/P_plen_nom;

M_res_nom   = -(M_N_used_nom*RN*Tc - P_plen_nom*V_plen_nom*k)/(2*RN*Tc);
P_tank_nom  = -(P_f_nom*(M_N_used_nom*RN*Tc + P_plen_nom*V_plen_nom*k))/(M_N_used_nom*RN*Tc - P_plen_nom*V_plen_nom*k);
V_tank_nom  = -(M_N_used_nom*RN*Tc - P_plen_nom*V_plen_nom*k)/(2*P_f_nom);
M_N_tot_nom = M_res_nom + M_N_used_nom;

%% Montecarlo


N = 10000;
rt_avg = sqrt(At/pi);
rt_std = 5e-6;
rt_val = normrnd(rt_avg, rt_std, [N,1]);
At_val = rt_val.^2.*pi;

r_inj_std = 5e-6;       
r_inj_avg = r_inj; 
r_inj_val = normrnd(r_inj_avg, r_inj_std, [N,1]);
A_inj_val = r_inj_val.^2.*pi;

K_tot_val = k_v*A_p/A_v + k_p + k_inj./A_inj_val*A_p;

re_avg = sqrt(Ae/pi);
re_std = 5e-6;
re_val = normrnd(re_avg, re_std, [N,1]);
Ae_val = re_val.^2.*pi;

toll = 1e-10;

T_avg = 0;
Isp_avg = 0;
P_c_mc_avg = 0;
m_dot_mc_avg = 0;
M_N_tot_avg = 0;
P_tank_avg = 0;
V_tank_avg = 0;
M_N_used_avg = 0;

for ii = 1:N
    
    A_tt  = At_val(ii);
    k_tot = K_tot_val(ii);
    Aee = Ae_val(ii);
    it = 0;

    P_cc = Pc;
    m_dot2 = 1;
    m_dot1 = 2;

    while(abs(m_dot2 - m_dot1) > toll && it < 1000)
        it = it + 1;

        m_dot1 =  P_cc*A_tt/c_star;
        v  = m_dot1 /( A_p* rho);
        deltaP= 0.5*rho*v^2*k_tot;
        P_cc = P_plen_nom - deltaP;
        
        m_dot2= P_cc*A_tt/c_star ;
       
        if it > 999
            disp('It did not converge')
        end   

    end

    eps = Aee/A_tt;
    
    aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
    bb = 2/k;
    cc = (k-1)/k;
    
    fun = @(x) x^bb*(1-x^cc)-aa;
    dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);
    
    c1 = 0;
    c2 = 0.5;
    tol = 1e-6;
    err = tol + 1;
    it = 0;
    
    while (it < 3 && err > tol ) 
        it=it+1;
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
    while (it< 5 && err> tol)
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
   
    Pee = xn * P_cc;

    P_c_mc(ii) = P_cc;
    m_dot_mc(ii) = m_dot2;
    m_dot_max = P_cc*A_tt*sqrt(k/RN/Tc)*(2/(k+1))^((k+1)/2/(k-1));
    m_dot_ratio(ii) = abs(m_dot_max-m_dot2)/m_dot_max;
    v_e = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pee/P_cc)^((k-1)/k)));
    T(ii) =  m_dot1*v_e + Aee*Pee;
    Isp_mc(ii) = T(ii)/(m_dot2*g);

    % mass
    MR_mc = exp(DeltaV/(Isp_mc(ii)*g));
    M_N_used_mc(ii) = Mass*(MR_mc-1);

    % pressure losses
    v_p_mc = m_dot_mc(ii)/(rho*A_p);
    Delta_P_tot_mc = 0.5*rho*v_p_mc^2*(k_p + A_p/A_v*k_v + A_p/A_inj*k_inj);
    Delta_P_plen_mc = 0.5*rho*v_p_mc^2*A_p/A_inj*k_inj;
    
    V_plen_mc = M_N_used_mc(end)*RN*Tc/P_plen_nom;
    
    M_res_mc(ii)   = -(M_N_used_mc(end)*RN*Tc - P_plen_nom*V_plen_mc*k)/(2*RN*Tc);
    P_tank_mc(ii)  = -(P_f_nom*(M_N_used_mc(end)*RN*Tc + P_plen_nom*V_plen_mc*k))/(M_N_used_mc(end)*RN*Tc - P_plen_nom*V_plen_mc*k);
    V_tank_mc(ii)  = -(M_N_used_mc(end)*RN*Tc - P_plen_nom*V_plen_mc*k)/(2*P_f_nom);
    M_N_tot_mc(ii) = M_res_mc(ii) + M_N_used_mc(end);


    T_avg(ii+1) = T_avg(ii) + 1/(ii)*(T(ii) - T_avg(ii));
    Isp_avg(ii+1) = Isp_avg(ii) + 1/(ii)*(Isp_mc(ii) - Isp_avg(ii));    
    P_c_mc_avg(ii+1) = P_c_mc(ii) + 1/(ii)*(P_c_mc(ii) - P_c_mc_avg(ii));
    m_dot_mc_avg(ii+1) = m_dot_mc_avg(ii) + 1/(ii)*(m_dot_mc(ii) - m_dot_mc_avg(ii));
    M_N_tot_avg(ii+1) = M_N_tot_avg(ii) + 1/(ii)*(M_N_tot_mc(ii) - M_N_tot_avg(ii));    
    
    M_N_used_avg(ii+1) = M_N_used_avg(ii) + 1/(ii)*(M_N_used_mc(end) - M_N_used_avg(ii));
    V_tank_avg(ii+1) = V_tank_avg(ii) + 1/(ii)*(V_tank_mc(ii) - V_tank_avg(ii));
    
    T_std(ii) = std(T(1:ii));
    Isp_std(ii) = std(Isp_mc(1:ii));

    P_c_mc_std(ii)   = std(P_c_mc(1:ii));
    m_dot_mc_std(ii) = std(m_dot_mc(1:ii));
    M_N_tot_std(ii)  = std(M_N_tot_mc(1:ii));

    M_N_used_std(ii) = std(M_N_used_mc(1:end)); 
    V_tank_std(ii)   = std(V_tank_mc(1:ii));
end

M_N_used_w = M_N_used_avg(end) + 3*M_N_used_std(end);

P_tank_w = (P_f_nom*(k + 1))/(k - 1);
V_tank_w = -(RN*Tc*(M_N_used_w - M_N_used_w*k))/(2*P_f_nom);
M_res_w  = (M_N_used_w*k)/2 - M_N_used_w/2;

figure 
subplot(2,2,1)
plot(T_avg(2:end))
%ylim([0 0.01])
title T_{avg}
subplot(2,2,2)
plot(T_std)
title T_{std}
subplot(2,2,3)
plot(Isp_avg(2:end))
title Isp_{avg}
subplot(2,2,4)
plot(Isp_std)
title Isp_{std}

figure 
subplot(2,2,1)
plot(M_N_used_avg(2:end))
title M_{N-used-avg}
subplot(2,2,2)
plot(M_N_used_std)
title M_{N-used-std}
subplot(2,2,3)
plot(V_tank_avg(2:end))
title V_{tank-avg}
subplot(2,2,4)
plot(V_tank_std)
title V_{tank-std}


%% leackage losses computed for 1 year
rho_stand = 1.176; %kg/m^3 cond standard
leack_rate = 10^-5; %scc/s standard cubic centieters every second
t_leack = 60*60*24*365;
tot_leack_scc = leack_rate*t_leack;
tot_mass_leack = tot_leack_scc*10^-6/rho_stand;
tot_N_mas = M_N_used_w + M_res_w - tot_mass_leack;
leacked_mass_ratio = tot_mass_leack/(M_N_used_w + M_res_w)*100;

% with a given volume calculate the loading the tank pressurization with a
% 20% margin

P_tank_marg = (M_N_used_w + M_res_w)*1.2*RN*Tc/V_tank_w
M_tot_marg  = P_tank_marg*V_tank_w/(RN*Tc) 
clc
clear
close all

%Constants

R = 8.314; 
g = 9.81;


%Initial conditions

T = 5e-3;           % Thrust (N)
Pc = 100000;        % chamber pressure (Pa) 1bar-4bar
Tc = 300;           % chamber temperature (K) fixed
MM = 0.028 ;        % kg/mol
k = 1.4;
eps = 4;            % area ratio from 4-20
DeltaV = 0.5;       % requested total DV (m/s)
Mass = 4;           % 4 U cubesat mass
rho = 807;

RN = R/MM;

% feeding system
k_p = 1.1;  % P drop coeff of the pipes
k_v = 1.3;  % P drop valves
k_inj = 0.6;    % P_drop injector

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

% Area Throat
At = (T/Pc)/(sqrt(2*(k^2/(k-1))*(2/(k+1))^((k+1)/(k-1)))*sqrt(1-(Pe/Pc)^((k-1)/k))+eps*(Pe-0)/Pc);
rt=sqrt(At/pi);
Dt=2*rt;

% Exit Velocity
v_exit = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pe/Pc)^((k-1)/k)));

% Exit Area
Ae = eps*At;
re=sqrt(Ae/pi);
De=2*re;

% Mass Flow Rate
m_dot=(T-Ae*(Pe-0))/v_exit;

% Isp
Isp_nom = T/m_dot/g;

%cstar
c_star = At*Pc/m_dot;

%% Preliminary sizing
% mass
MR = exp(DeltaV/(Isp_nom*g));
M_N_used = Mass*(MR-1);

% pressure losses
v_p = m_dot/(rho*A_p);
Delta_P_tot = 0.5*rho*v_p^2*(k_p + A_p/A_v*k_v + A_p/A_inj*k_inj);
Delta_P_plen = 0.5*rho*v_p^2*A_p/A_inj*k_inj;
P_plen = Pc + Delta_P_plen;
P_f = Pc + Delta_P_tot;
P_f = P_f*1.2;

% imposed volume of tank
%V_tank = 1e-4; % thesis Cold_gas Marcus lore
%M_N_residual = V_tank*P_f/(R*Tc); 
%P_tank = (M_N_residual+M_N_used)*R*Tc/V_tank;

V_plen = M_N_used*RN*Tc/P_plen;

eq =@(V_tank) V_tank.*P_f./RN./Tc + M_N_used - P_f.*P_plen.*V_plen./RN./Tc.*k./(2.*P_f+M_N_used.*RN*Tc./V_tank);

V_tank = fsolve(eq,1e-4);
P_tank = (M_N_used + P_f*V_tank/RN/Tc)*RN*Tc/V_tank;

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

toll = 1e-8;

T = [];
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
        P_cc = P_plen - deltaP;
        
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

    v_e = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pee/P_cc)^((k-1)/k)));
    T(ii) =  m_dot1*v_e + Aee*Pee;
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
ylim([0 0.01])
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
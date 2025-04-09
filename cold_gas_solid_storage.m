clc
clear
close all

%Constants

R = 8.314; 
g = 9.81;


%Initial conditions

T_nom = 5e-3;           % Thrust (N)
% Pressure in the plenum mid, (1 - 4.5)
%P_pl_mid = 260000;  % for complete discharge
P_pl_mid = 235600;
Tc = 300;           % chamber temperature (K) fixed
MM = 0.028 ;        % kg/mol
k = 1.4;
eps = 4;            % area ratio from 4-20
DeltaV = 0.5;       % requested total DV (m/s)
Mass = 4;           % 4 U cubesat mass
rho = 807;
t_burn = DeltaV/(T_nom/Mass);
P_tank_nom = 460000;


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

v_p = 1.2;   % guessed 

Delta_P_plen_nom = 0.5*rho*v_p^2*k_inj;
Pc = P_pl_mid - Delta_P_plen_nom;

%% Nominal Design


aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
bb = 2/k;
cc = (k-1)/k;

fun = @(x) x^bb*(1-x^cc)-aa;
dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);

c1 = 0;
c2 = 0.5;
tol = 1e-9;
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
At = (T_nom/Pc)/(sqrt(2*(k^2/(k-1))*(2/(k+1))^((k+1)/(k-1)))*sqrt(1-(Pe/Pc)^((k-1)/k))+eps*(Pe-0)/Pc);
rt=sqrt(At/pi);
Dt=2*rt;

% Exit Velocity
v_exit = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pe/Pc)^((k-1)/k)));


% Exit Area
Ae = eps*At;
re=sqrt(Ae/pi);
De=2*re;

% Mass Flow Rate
m_dot_nom=(T_nom-Ae*(Pe-0))/v_exit;

% Isp
Isp_nom = T_nom/m_dot_nom/g;

%cstar
c_star = At*Pc/m_dot_nom;


% variable fixed all noozzle cstar T 
% we bring as guesses Pc mid 

charge_mass = 0.0003; % preliminar guess
Delta_p_mid = 360000;
V_plen = charge_mass*RN*Tc/Delta_p_mid;


step = 0.1;
P_tank = P_tank_nom;
ii = 1;
toll = 1e-8;
M_used = 0;
while P_tank(end) > 100000
    m_dot2 = 5;
    m_dot1 = 19;
    it = 0;
    P_cc = Pc;
    while(abs(m_dot2 - m_dot1) > toll && it < 10000)
        it = it + 1;

        m_dot1 =  P_cc*At/c_star;
        v  = m_dot1 /( A_inj* rho);
        deltaP= 0.5*rho*v^2*k_inj;
        P_cc = P_tank(end) - deltaP;
        
        m_dot2= P_cc*At/c_star ;

        if it > 9999
            disp('Il calcolo non converge')
        end   
    end

    deltaP_tank = m_dot1*step*RN*Tc/V_plen;
    P_tank = [P_tank , P_tank(end)-deltaP_tank];

    m_dot_vect(ii) = m_dot2;
    
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
    Pee = xn * P_cc;

    v_e     = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pee/P_cc)^((k-1)/k)));
    T(ii)   = m_dot1*v_e + Ae*Pee;
    Isp(ii) = T(ii)/(m_dot2*g);
    M_used  = M_used + m_dot2*step;
    ii = ii + 1;
end

T_avg_nom  = mean(T)
t_burn_each = length(T)*step

figure
plot(T)
%% MonteCarlo

N = 10000;

rt_avg = sqrt(At/pi);
rt_std = 5e-6;
rt_val = normrnd(rt_avg, rt_std, [N,1]);
At_val = rt_val.^2.*pi;

r_inj_std = 5e-6;       
r_inj_avg = r_inj; 
r_inj_val = normrnd(r_inj_avg, r_inj_std, [N,1]);
A_inj_val = r_inj_val.^2.*pi;

re_avg = sqrt(Ae/pi);
re_std = 5e-6;
re_val = normrnd(re_avg, re_std, [N,1]);
Ae_val = re_val.^2.*pi;

P_tank_val = normrnd(P_tank_nom, P_tank_nom*0.01, [N,1]);

toll = 1e-8;

T_avg = 0;
t_b_each_avg = 0;
for jj = 1:N

    A_tt = At_val(jj);
    Aee = Ae_val(jj);
    A_injj = A_inj_val(jj);
    it = 0;

    step = 0.1;
    P_tank = P_tank_val(jj);
    T = [];
    ii = 1;
    toll = 1e-8;
    while P_tank(end) > 100000
        m_dot2 = 5;
        m_dot1 = 19;
        it = 0;
        P_cc = Pc;
        while(abs(m_dot2 - m_dot1) > toll && it < 10000)
            it = it + 1;
    
            m_dot1 =  P_cc*A_tt/c_star;
            v  = m_dot1 /( A_injj* rho);
            deltaP= 0.5*rho*v^2*k_inj;
            P_cc = P_tank(end) - deltaP;
            
            m_dot2= P_cc*A_tt/c_star ;
    
            if it > 9999
                disp('Il calcolo non converge')
            end   
        end
    
        deltaP_tank = m_dot1*step*RN*Tc/V_plen;
        P_tank = [P_tank , P_tank(end)-deltaP_tank];
    
        m_dot_vect(ii) = m_dot2;
        
        aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
        bb = 2/k;
        cc = (k-1)/k;
        
        fun = @(x) x^bb*(1-x^cc)-aa;
        dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);
        
        c1 = 0;
        c2 = 0.5;
        tol = 1e-9;
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
        while (it< 20 && err> tol)
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
        Pee = xn * P_cc;
    
        v_e     = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pee/P_cc)^((k-1)/k)));
        T(ii)   = m_dot1*v_e + Ae*Pee;
        Isp(ii) = T(ii)/(m_dot2*g);
        ii = ii + 1;
    end
    T_vect(jj)          = mean(T);
    t_b_each(jj)        = step*length(T);
    T_avg(jj+1)         = T_avg(jj) + 1/jj*(T_vect(jj) - T_avg(jj));
    t_b_each_avg(jj+1) = t_b_each_avg(jj) + 1/jj*(t_b_each(jj) - t_b_each_avg(jj));
    T_std(jj) = std(T_vect(1:jj));
end

%% plots

figure
subplot(2,2,1)
plot(T_avg(2:end))
subplot(2,2,2)
plot(T_std(2:end))
subplot(2,2,3)
plot(t_b_each_avg(2:end))
title t_{b-each-avg}
%% edge case
% the first edge case investigated is the first charge, that does not reach
% the nominal high pressure

P_plen_first_charge = charge_mass*RN*Tc/V_plen; % sono un idiota è quello prima meno 1

P_tank_fc = P_plen_first_charge;
ii = 1;
toll = 1e-9;

while P_tank_fc(end) > 100000
    m_dot2 = 5;
    m_dot1 = 19;
    it = 0;
    P_cc = Pc;
    while(abs(m_dot2 - m_dot1) > toll && it < 10000)
        it = it + 1;

        m_dot1 =  P_cc*At/c_star;
        v  = m_dot1 /( A_inj* rho);
        deltaP= 0.5*rho*v^2*k_inj;
        P_cc = P_tank_fc(end) - deltaP;
        
        m_dot2= P_cc*At/c_star ;

        if it > 9999
            disp('Il calcolo non converge')
        end   
    end

    deltaP_tank = m_dot1*step*RN*Tc/V_plen;
    P_tank_fc = [P_tank_fc , P_tank_fc(end)-deltaP_tank];

    m_dot_vect_fc(ii) = m_dot2;
    
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
    Pee = xn * P_cc;

    v_e        = sqrt(2*k/(k-1)*R/MM*Tc*(1-(Pee/P_cc)^((k-1)/k)));
    T_fc(ii)   = m_dot1*v_e + Ae*Pee;
    Isp_fc(ii) = T_fc(ii)/(m_dot2*g);

    ii = ii + 1;
end

T_avg_fc  = mean(T_fc);
t_burn_fc = length(T_fc)*step;


%% N charges

Delta_V_fc = T_avg_fc/Mass*t_burn_fc
Delta_V_each = T_avg(end)/Mass*t_b_each_avg(end)

N_charges = (DeltaV - Delta_V_fc)/Delta_V_each + 1

mass_totale = N_charges*charge_mass

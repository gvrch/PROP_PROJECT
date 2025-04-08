clear 
close all
clc
%% NOMINAL

% DATA
eps   = 9;            % Rapporto tra le aree (ugello convergente-divergente)  [adim]
k       = 1.327;        % Rapporto dei calori specifici                         [adim]
g       = 9.80665;         % Gravity acceleration                                  [m/s^2]
Pc_d    = 230000;       % Evaporation chamber Pressure                          [Pa]
R       = 8.314;        % Gas Constant                                          [J/mol*K]
Mmol    = 0.018;        % Water Molar Mass                                      [kg/mol]
rho_l   = 999;          % Density liquid water                                  [kg/m^3]
cp_l    = 4.1838e3;     % Specific Heat liquid water                            [J/kgK]
cp_v    = 2.0256e3;     % Specific Heat gas water                               [J/kgK]
lambda  = 2260e+3;      % Water latent evaporation heat                         [J/kg]
T_eb    = 124 + 273.15; % Temperature of ebollition                             [K]
T_inj   = 25 + 273.15;  % temperature of Injection                              [K]
alpha   = deg2rad(30);  %                                                       [deg]
beta    = deg2rad(45);  %                                                       [deg]
Cd      = 0.7;          % Drag coefficient                                      [adim]
dP_inj_perc = 0.1;      % Pressure lost for injection   !!!CHECK!!!
T     = 5e-3;         % Wanted Thrust                                         [N]
dv      = 0.5;          % Wanted Delta v                                        [m/s]
mass    = 4;            % S/C mass                                              [kg]
N = 5000;              % N MONTECARLO
% GAMMA = sqrt(k*(2/(k+1))^((k+1)/(k-1)));
% Injection - Cd (TBR)
k_inj = 1/Cd^2;

% Calculate the Pressure Ratio (Pe/Pc)
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
p_ratio = xn;
Pe_d = p_ratio*Pc_d;

% Compute throat and exit area
CT = k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-p_ratio^((k-1)/k)))+p_ratio*eps;
At = T/Pc_d/CT;
r_t = sqrt(At/pi);
Ae = eps*At;
r_e = sqrt(Ae/pi);

% Dimension nozzle
eps_conv = 1/0.2 * sqrt((2/(k+1)*(1+(k-1)/2*0.2^2))^((k+1)/(k-1)));
Ac = eps_conv * At;
r_c = sqrt(Ac/pi);

L_conv = 1/2*(2*(r_c-r_t)/tan(beta)); % ipotesi che Ae è uagule a quwlla di ingresso nel convergente
L_div = 1/2*(2*(r_e-r_t)/tan(alpha));
L_nozzle = L_conv + L_div;
l = 1/2*(1+cos(alpha));
l = 1;

% Define some usefull constant value 
K       = sqrt(2*k/(k-1) * R/Mmol * (1 - p_ratio^((k-1)/k))); % V_e = k * sqrt(Tc)
aa      = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v;
cc      = cp_v*((T-Pe_d*Ae)/K)^2;
Q_min   = sqrt(4*aa*cc);

Q_design = linspace(Q_min, Q_min+1.5, 100);   % Only for design

unstart_design = zeros(length(Q_design),1);
choke_design = zeros(length(Q_design),1);
T_avg_design = zeros(length(Q_design),1);
T_std_design = zeros(length(Q_design),1);
Pc_avg_design = zeros(length(Q_design),1);
Pc_std_design = zeros(length(Q_design),1);
Pe_avg_design = zeros(length(Q_design),1);
Pe_std_design = zeros(length(Q_design),1);
m_dot_avg_design = zeros(length(Q_design),1);
m_dot_std_design = zeros(length(Q_design),1);
ve_avg_design = zeros(length(Q_design),1);
ve_std_design = zeros(length(Q_design),1);
Tc_avg_design = zeros(length(Q_design),1);
Tc_std_design = zeros(length(Q_design),1);
Isp_avg_design = zeros(length(Q_design),1);
Isp_std_design = zeros(length(Q_design),1);

for iii = 1:length(Q_design)
    Q_dot = Q_design(iii);
    
    % Study the flow in function of Q_dot
    
    % Define some usefull constant value 
    K       = sqrt(2*k/(k-1) * R/Mmol * (1 - p_ratio^((k-1)/k))); % V_e = k * sqrt(Tc)
    aa      = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v;
    cc      = cp_v*((T-Pe_d*Ae)/K)^2;
    Q_min   = sqrt(4*aa*cc);
        
    % mass flow rate
   % m_dot = (Q_dot + sqrt(Q_dot^2 - Q_min^2))/(2*aa);  % LOW TEMP BRANCH
     m_dot = (Q_dot - sqrt(Q_dot^2 - Q_min^2))/(2*aa);  % HIGH TEMP BRANCH
    
    % Combustion Chamber Temperature
    Tc = ((T-Pe_d*Ae)/(m_dot *K))^2;
    
    % Exit Velocity
    v_e = K *sqrt(Tc);
    
    % Exit Temperature
    T_e = Tc*p_ratio^((k-1)/k);
    
    % Specific Impulse
    Isp = T/(m_dot*g);
    
    % Injection Area
    dP_inj = dP_inj_perc*Pc_d;
    Ainj = m_dot/(Cd*sqrt(2*dP_inj*rho_l));
    r_inj = sqrt(Ainj/pi);

    % dP_inj= 0.5*k_inj*m_dot^2/Ainj^2/rho_l;
    
    P_plenum = Pc_d + dP_inj;
    
    K_tot = 0.5/rho_l*(k_inj/Ainj^2);
    dP_tot = K_tot*m_dot^2;
    
    %% MONTECARLO
    
    % Definizione parametri
    tol = 1e-12;         % Tolleranza convergenza
    
    r_inj_std = (0.5e-6)/3; % TOLLERANZA A 3 sigma (1mm su diametro) ----------------- CITARE ANDREA
    r_inj_avg = r_inj;
    r_t_std = (0.5e-6)/3;
    r_t_avg = r_t;
    Q_dot_std = Q_dot*2/100/3;
    Q_dot_avg = Q_dot;
    
    % Ainj_std = 0;
    % At_std = 0;
    % Q_dot_std = 0;
    
    r_inj_val = normrnd(r_inj_avg, r_inj_std, [N,1]);
    r_t_val = normrnd(r_t_avg, r_t_std, [N,1]);
    Q_dot_val =  normrnd(Q_dot_avg, Q_dot_std, [N, 1]);
    Ainj_val = pi*r_inj_val.^2;
    At_val = pi*r_t_val.^2;

    K_tot_val = 0.5/rho_l.*(k_inj./Ainj_val.^2);
    % dP_tot = K_tot*m_dot^2;
        
    T_val = zeros(N,1);
    Pc_val = zeros(N,1);
    Pe_val = zeros(N,1);
    m_dot_val = zeros(N,1);
    ve_val = zeros(N,1);
    Tc_val = zeros(N,1);
    Isp_val = zeros(N,1);
    
    T_avg = zeros(N+1,1);
    T_std = zeros(N+1,1);
    Pc_avg = zeros(N+1,1);
    Pc_std = zeros(N+1,1);
    Pe_avg = zeros(N+1,1);
    Pe_std = zeros(N+1,1);
    m_dot_avg = zeros(N+1,1);
    m_dot_std = zeros(N+1,1);
    ve_avg = zeros(N+1,1);
    ve_std = zeros(N+1,1);
    Tc_avg = zeros(N+1,1);
    Tc_std = zeros(N+1,1);
    Isp_avg = zeros(N+1,1);
    Isp_std = zeros(N+1,1);
    
    
    unstart = 0;
    choke = 0;
    for i = 1:N
        Ainj  = Ainj_val(i);
        At    = At_val(i);
        K_tot = K_tot_val(i); % dP_tot = K_tot*m_dot^2;
        Q_dot = Q_dot_val(i);
    
        eps = Ae/At;
    
        % Calculate the Pressure Ratio (Pe/Pc)
        aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
        bb = 2/k;
        cc = (k-1)/k;
        
        % Resolve with Newton method
        fun  = @(x) x^bb*(1-x^cc)-aa;
        dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);
        
        c1  = 0;
        c2  = 0.5;
        tol = 1e-16;
        err = tol + 1;
        it  = 0;
        
        while (it < 5 && err > tol ) 
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
        while (it < 200 && err> tol)
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
        % Pe = xn * Pc;
        p_ratio_i = xv;
        % K_ratio2= k*(2/(k+1))^((k+1)/(k-1))*Mmol/R;
        M =0.05;
        GAMMA = 1/(sqrt(k)*M*(1+(k-1)/2*M^2)^((1+k)/(2-2*k)));
        K_ratio2= GAMMA^2*Mmol/R;
    
        AA = cp_v*K_ratio2*At^2*K_tot^2;
        BB = 0;
        CC = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v -2*P_plenum*K_tot*cp_v*K_ratio2*At^2;
        DD = - Q_dot;
        EE = cp_v*K_ratio2*P_plenum^2*At^2;

        Delta = 256*AA^3*EE^3 - 128*AA^2*CC^2*EE^2 + 144*AA^2*CC*DD^2*EE - 27*AA^2*DD^4 + ...
            16*AA*CC^4*EE - 4*AA*CC^3*DD^2;
        P = 8*AA*CC;
        D = 64*AA^3*EE - 16*AA^2*CC^2;

        Delta0 = CC^2 + 12*AA*EE;
        Delta1 = 2*CC^3 + 27*BB^2*EE + 27*AA*DD^2 - 72*AA*CC*EE;

        q = CC/AA;
        r = DD/AA;
        s = EE/AA;ij

        if Delta > 0
            sucaaaa = 0;
        end

        fun = @(x) AA*x.^4 + CC*x.^2 + DD*x + EE;
        dfun = @(x) 4*AA*x.^3 + 2*CC*x + DD;
        dfun2 = @(x) 12*AA*x.^2 + 2*CC;
    
        toll = 1e-12;
        err = toll + 1;
        iter = 0;
        xv = m_dot;
        while (iter < 20 && err> toll)
           dfx = dfun2(xv);
           if dfx == 0
              error(' Arresto per azzeramento di dfun');
           else
              xn = xv - dfun(xv)/dfx;
              err = abs(dfun(xn));
              iter = iter+1;
              xv = xn;
           end
        end
        xmin = xv;
    
        c1 = xmin;      % LOW TEMP BRANCH
        c2 = 1e-4;
        % c1 = 1e-7;    % HIGH TEMP BRANCH
        % c2 = xmin;
        toll = 1e-10;
        err = toll + 1;
        iter = 0;
        while (iter < 40 && err > toll ) 
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
        m_dot = x;
    
        Pc = P_plenum- K_tot*m_dot^2;
        Tc = (Q_dot-m_dot*(cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v))/m_dot/cp_v;
        unstart = unstart + (Tc < T_eb);
        Pe = Pc*p_ratio_i;
        ve = sqrt(2*k/(k-1)*R/Mmol*Tc*(1-(Pe/Pc)^((k-1)/k)));
    
        m_dot_max = Pc*At*sqrt(k/(R/Mmol)/Tc)*(2/(k+1))^((k+1)/2/(k-1));
        choke = choke + ((m_dot-m_dot_max)/m_dot > 1e-3);


        c_star1 = Pc*At/m_dot;
        c_star2 = sqrt(R/Mmol*Tc)/GAMMA;
        % choke = choke + (abs(c_star1-c_star2)/c_star2 > 1e-4);


        T_val(i)=(m_dot*l*ve+Pe*Ae)*(Tc > T_eb);
        Pc_val(i)=Pc*(Tc > T_eb);
        Pe_val(i)=Pe*(Tc > T_eb);
        m_dot_val(i)=m_dot*(Tc > T_eb);
        ve_val(i)=ve*(Tc > T_eb);
        Isp_val(i) = T_val(i)/m_dot_val(i)/g*(Tc > T_eb);
        Tc_val(i) = Tc*(Tc > T_eb);
    
        T_avg(i+1)=(T_avg(i)*(i-1)+T_val(i))/i;
        T_std(i+1)=std(T_val(1:i));
        Pc_avg(i+1)=(Pc_avg(i)*(i-1)+Pc_val(i))/i;
        Pc_std(i+1)=std(Pc_val(1:i));
        Pe_avg(i+1)=(Pe_avg(i)*(i-1)+Pe_val(i))/i;
        Pe_std(i+1)=std(Pe_val(1:i));
        m_dot_avg(i+1)=(m_dot_avg(i)*(i-1)+m_dot_val(i))/i;
        m_dot_std(i+1)=std(m_dot_val(1:i));
        ve_avg(i+1)=(ve_avg(i)*(i-1)+ve_val(i))/i;
        ve_std(i+1)=std(ve_val(1:i));
        Tc_avg(i+1)=(Tc_avg(i)*(i-1)+Tc_val(i))/i;
        Tc_std(i+1)=std(Tc_val(1:i));
        Isp_avg(i+1)=(Isp_avg(i)*(i-1)+Isp_val(i))/i;
        Isp_std(i+1)=std(Isp_val(1:i));

    end
    
    start_frac = 1 - unstart/N;
    start_frac = start_frac + (start_frac < 1e-3);
    unstart_design(iii) = unstart/N;
    choke_design(iii) = choke/N;
    T_avg_design(iii) = T_avg(end)/start_frac;
    T_std_design(iii) = T_std(end)/start_frac;
    Pc_avg_design(iii) = Pc_avg(end)/start_frac;
    Pc_std_design(iii) = Pc_std(end)/start_frac;
    Pe_avg_design(iii) = Pe_avg(end)/start_frac;
    Pe_std_design(iii) = Pe_std(end)/start_frac;
    m_dot_avg_design(iii) = m_dot_avg(end)/start_frac;
    m_dot_std_design(iii) = m_dot_std(end)/start_frac;
    ve_avg_design(iii) = ve_avg(end)/start_frac;
    ve_std_design(iii) = ve_std(end)/start_frac;
    Tc_avg_design(iii) = Tc_avg(end)/start_frac;
    Tc_std_design(iii) = Tc_std(end)/start_frac;
    Isp_avg_design(iii) = Isp_avg(end)/start_frac;
    Isp_std_design(iii) = Isp_std(end)/start_frac;
end

%% PLOTS
figure(1)
plot(Q_design,unstart_design*100,Q_design,choke_design*100)
legend("unstart","choke")

figure(2)
subplot(2,7,1)
plot(Q_design,T_avg_design)
legend("T_avg")
subplot(2,7,8)
plot(Q_design,T_std_design)
legend("T_std")
subplot(2,7,2)
plot(Q_design,Pc_avg_design)
legend("Pc_avg")
subplot(2,7,9)
plot(Q_design,Pc_std_design)
legend("Pc_std")
subplot(2,7,3)
plot(Q_design,Pe_avg_design)
legend("Pe_avg")
subplot(2,7,10)
plot(Q_design,Pe_std_design)
legend("Pe_std")
subplot(2,7,4)
plot(Q_design,m_dot_avg_design)
legend("m_dot_avg")
subplot(2,7,11)
plot(Q_design,m_dot_std_design)
legend("m_dot_std")
subplot(2,7,5)
plot(Q_design,ve_avg_design)
legend("ve_avg")
subplot(2,7,12)
plot(Q_design,ve_std_design)
legend("ve_std")
subplot(2,7,6)
plot(Q_design,Tc_avg_design)
legend("Tc_avg")
subplot(2,7,13)
plot(Q_design,Tc_std_design)
legend("Tc_std")
subplot(2,7,7)
plot(Q_design,Isp_avg_design)
legend("Isp_avg")
subplot(2,7,14)
plot(Q_design,Isp_std_design)
legend("Isp_std")

%% Most Reliable
[unstart_best, iii] = min(unstart_design);

Q_dot = Q_design(iii);

% Study the flow in function of Q_dot

% Define some usefull constant value 
K       = sqrt(2*k/(k-1) * R/Mmol * (1 - p_ratio^((k-1)/k))); % V_e = k * sqrt(Tc)
aa      = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v;
cc      = cp_v*((T-Pe_d*Ae)/K)^2;
Q_min   = sqrt(4*aa*cc);
    
% mass flow rate
m_dot = (Q_dot + sqrt(Q_dot^2 - Q_min^2))/(2*aa)  % LOW TEMP BRANCH
% m_dot = (Q_dot - sqrt(Q_dot^2 - Q_min^2))/(2*aa)  % HIGH TEMP BRANCH

% Combustion Chamber Temperature
Tc = ((T-Pe_d*Ae)/(m_dot *K))^2

% Exit Velocity
v_e = K *sqrt(Tc)

% Exit Temperature
T_e = Tc*p_ratio^((k-1)/k)

% Specific Impulse
Isp = T/(m_dot*g)

% Injection Area
dP_inj = dP_inj_perc*Pc_d
Ainj = m_dot/(Cd*sqrt(2*dP_inj*rho_l));
r_inj = sqrt(Ainj/pi);

% dP_inj= 0.5*k_inj*m_dot^2/Ainj^2/rho_l;

P_plenum = Pc_d + dP_inj

K_tot = 0.5/rho_l*(k_inj/Ainj^2);
dP_tot = K_tot*m_dot^2;

% MONTECARLO

% Definizione parametri
tol = 1e-12;         % Tolleranza convergenza

r_inj_std = (0.5e-6)/3; % TOLLERANZA A 3 sigma (1mm su diametro) ----------------- CITARE ANDREA
r_inj_avg = r_inj;
r_t_std = (0.5e-6)/3;
r_t_avg = r_t;
Q_dot_std = Q_dot*2/100/3;
Q_dot_avg = Q_dot;

% Ainj_std = 0;
% At_std = 0;
% Q_dot_std = 0;

r_inj_val = normrnd(r_inj_avg, r_inj_std, [N,1]);
r_t_val = normrnd(r_t_avg, r_t_std, [N,1]);
Q_dot_val =  normrnd(Q_dot_avg, Q_dot_std, [N, 1]);
Ainj_val = pi*r_inj_val.^2;
At_val = pi*r_t_val.^2;

K_tot_val = 0.5/rho_l.*(k_inj./Ainj_val.^2);
% dP_tot = K_tot*m_dot^2;
    
T_val = zeros(N,1);
Pc_val = zeros(N,1);
Pe_val = zeros(N,1);
m_dot_val = zeros(N,1);
ve_val = zeros(N,1);
Tc_val = zeros(N,1);
Isp_val = zeros(N,1);

T_avg = zeros(N+1,1);
T_std = zeros(N+1,1);
Pc_avg = zeros(N+1,1);
Pc_std = zeros(N+1,1);
Pe_avg = zeros(N+1,1);
Pe_std = zeros(N+1,1);
m_dot_avg = zeros(N+1,1);
m_dot_std = zeros(N+1,1);
ve_avg = zeros(N+1,1);
ve_std = zeros(N+1,1);
Tc_avg = zeros(N+1,1);
Tc_std = zeros(N+1,1);
Isp_avg = zeros(N+1,1);
Isp_std = zeros(N+1,1);

unstart = 0;
for i = 1:N
    Ainj  = Ainj_val(i);
    At    = At_val(i);
    K_tot = K_tot_val(i); % dP_tot = K_tot*m_dot^2;
    Q_dot = Q_dot_val(i);

    eps = Ae/At;

    % Calculate the Pressure Ratio (Pe/Pc)
    aa = 1/eps^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
    bb = 2/k;
    cc = (k-1)/k;
    
    % Resolve with Newton method
    fun  = @(x) x^bb*(1-x^cc)-aa;
    dfun = @(x) -x^(bb-1)*((cc+bb)*x^cc-bb);
    
    c1  = 0;
    c2  = 0.5;
    tol = 1e-7;
    err = tol + 1;
    it  = 0;
    
    while (it < 5 && err > tol ) 
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
    while (it < 15 && err> tol)
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
    % Pe = xn * Pc;
    p_ratio_i = xv;
    K_ratio2= k*(2/(k+1))^((k+1)/(k-1))*Mmol/R;

    AA = cp_v*K_ratio2*At^2*K_tot^2;
    BB = 0;
    CC = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v -2*P_plenum*K_tot*cp_v*K_ratio2*At^2;
    DD = - Q_dot;
    EE = cp_v*K_ratio2*P_plenum^2*At^2;

    fun = @(x) AA*x.^4 + CC*x.^2 + DD*x + EE;
    dfun = @(x) 4*AA*x.^3 + 2*CC*x + DD;
    dfun2 = @(x) 12*AA*x.^2 + 2*CC;

    toll = 1e-12;
    err = toll + 1;
    iter = 0;
    xv = m_dot;
    while (iter < 20 && err> toll)
       dfx = dfun2(xv);
       if dfx == 0
          error(' Arresto per azzeramento di dfun');
       else
          xn = xv - dfun(xv)/dfx;
          err = abs(dfun(xn));
          iter = iter+1;
          xv = xn;
       end
    end
    xmin = xv;

    c1 = xmin;      % LOW TEMP BRANCH
    c2 = 1e-4;
    % c1 = 1e-7;    % HIGH TEMP BRANCH
    % c2 = xmin;
    toll = 1e-10;
    err = toll + 1;
    iter = 0;
    while (iter < 30 && err > toll ) 
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
    m_dot = x;

    Pc = P_plenum- K_tot*m_dot^2;
    Tc = (Q_dot-m_dot*(cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v))/m_dot/cp_v;

    unstart = unstart + (Tc < T_eb);
    Tc = Tc*(Tc > T_eb);
    Pe = Pc*p_ratio_i;
    ve = sqrt(2*k/(k-1)*R/Mmol*Tc*(1-p_ratio_i^((k-1)/k)));

    if Tc > T_eb
        index = i - unstart;

        T_val(index)=(m_dot*l*ve+Pe*Ae)*(Tc > T_eb);
        Pc_val(index)=Pc*(Tc > T_eb);
        Pe_val(index)=Pe*(Tc > T_eb);
        m_dot_val(index)=m_dot*(Tc > T_eb);
        ve_val(index)=ve*(Tc > T_eb);
        Isp_val(index) = T_val(index)/m_dot_val(index)/g;
        Tc_val(index) = Tc;

        T_avg(index+1)=(T_avg(index)*(index-1)+T_val(index))/index;
        T_std(index+1)=std(T_val(1:index));
        Pc_avg(index+1)=(Pc_avg(index)*(index-1)+Pc_val(index))/index;
        Pc_std(index+1)=std(Pc_val(1:index));
        Pe_avg(index+1)=(Pe_avg(index)*(index-1)+Pe_val(index))/index;
        Pe_std(index+1)=std(Pe_val(1:index));
        m_dot_avg(index+1)=(m_dot_avg(index)*(index-1)+m_dot_val(index))/index;
        m_dot_std(index+1)=std(m_dot_val(1:index));
        ve_avg(index+1)=(ve_avg(index)*(index-1)+ve_val(index))/index;
        ve_std(index+1)=std(ve_val(1:index));
        Tc_avg(index+1)=(Tc_avg(index)*(index-1)+Tc_val(index))/index;
        Tc_std(index+1)=std(Tc_val(1:index));
        Isp_avg(index+1)=(Isp_avg(index)*(index-1)+Isp_val(index))/index;
        Isp_std(index+1)=std(Isp_val(1:index));
    end
end


%% Best Plots

figure(3)
subplot(2,7,1)
plot(T_avg)
legend("T_avg")
subplot(2,7,8)
plot(T_std)
legend("T_std")
subplot(2,7,2)
plot(Pc_avg)
legend("Pc_avg")
subplot(2,7,9)
plot(Pc_std)
legend("Pc_std")
subplot(2,7,3)
plot(Pe_avg)
legend("Pe_avg")
subplot(2,7,10)
plot(Pe_std)
legend("Pe_std")
subplot(2,7,4)
plot(m_dot_avg)
legend("m_dot_avg")
subplot(2,7,11)
plot(m_dot_std)
legend("m_dot_std")
subplot(2,7,5)
plot(ve_avg)
legend("ve_avg")
subplot(2,7,12)
plot(ve_std)
legend("ve_std")
subplot(2,7,6)
plot(Tc_avg)
legend("Tc_avg")
subplot(2,7,13)
plot(Tc_std)
legend("Tc_std")
subplot(2,7,7)
plot(Isp_avg)
legend("Isp_avg")
subplot(2,7,14)
plot(Isp_std)
legend("Isp_std")
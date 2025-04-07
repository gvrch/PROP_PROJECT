clear 
close all
clc
%% NOMINAL

addpath utils\
% DATA
eps_d     = 9;            % Rapporto tra le aree (ugello convergente-divergente)  [adim]
k       = 1.327;        % Rapporto dei calori specifici                         [adim]
g       = 9.80665;         % Gravity acceleration                                  [m/s^2]
Pc_d      = 230000;       % Evaporation chamber Pressure                          [Pa]
R       = 8.314;        % Gas Constant                                          [J/mol*K]
Mmol    = 0.018;        % Water Molar Mass                                      [kg/mol]
rho_l   = 1000;         % Density liquid water                                  [kg/m^3]
rho_v   = 0.6;          % Density gas water                                     [kg/m^3]
cp_l    = 4.1838e3;     % Specific Heat liquid water                            [J/kgK]
cp_v    = 2.0256e3;     % Specific Heat gas water                               [J/kgK]
lambda  = 2260e+3;      % Water latent evaporation heat                         [J/kg]
T_eb    = 124 + 273.15;    % Temperature of ebollition                             [K]
T_inj   = 25 + 273.15;  % temperature of Injection                              [K]
alpha   = deg2rad(30);  %                                                       [deg]
beta    = deg2rad(45);  %                                                       [deg]
Cd      = 0.7;          % Drag coefficient                                      [adim]
dP_inj_perc = 0.1;      % Pressure lost for injection   !!!CHECK!!!
T_d       = 5e-3;         % Wanted Thrust                                         [N]
dv      = 0.5;          % Wanted Delta v                                        [m/s]
mass    = 4;            % S/C mass                                              [kg]
Q_design = linspace(11.8, 13, 20);           %   SCELTO DAL NOMINALE PLOT
N = 10000;              % N MONTECARLO

k_pipe = 1;        % Coefficiente di perdita tubi
k_v = 1.3;      % Coefficiente di perdita valvola
% k_inj = 1.17;     % Coefficiente di perdita iniettore
k_inj = 1/Cd^2;

r_pipe = 500e-6;      % Area del tubo di alimentazione (m^2) CAMBIARE
Av = 32e-5;        % Area delle valvole (m^2)
Apipe = r_pipe^2*pi;

% Scraping NIST data
% load("H2O_data.mat")
% P = linspace(0.05, 1, 100);
% % P = 0.05;
% T = linspace(275, 1450, 100);
% f = @(x1, x2) interplinear(H2O.enthalpy, [x1 x2]);
% F = zeros(max(size(T)), max(size(P)));
% 
% for ii = 1:length(T)
%     if T(ii)>1445
%         T(ii)
%     end
%     for jj = 1:length(P)
%        F(ii, jj) = f(T(ii), P(jj));
%     end
% end
% 
% figure
% surf(F)

unstart_design = zeros(length(Q_design),1);
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
    eps = eps_d;
    Pc = Pc_d;
    T = T_d;
    
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
    
    % Compute exit pressure value
    Pe = xn * Pc;
    
    % Compute throat and exit area 
    
    CT = k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k)))+Pe/Pc*eps;
    At = T/Pc/CT;
    r_t = sqrt(At/pi);
    Ae = eps*At;
    r_e = sqrt(Ae/pi);
    
    % Dimensionate nozzle
    eps_conv = 1/0.2 * sqrt((2/(k+1)*(1+(k-1)/2*0.2^2))^((k+1)/(k-1)));
    Ac = eps_conv * At;
    r_c = sqrt(Ac/pi);
    
    L_conv = 1/2*(2*(r_c-r_t)/tan(beta)); % ipotesi che Ae è uagule a quwlla di ingresso nel convergente
    L_div = 1/2*(2*(r_e-r_t)/tan(alpha));
    L_nozzle = L_conv + L_div;
    l = 1/2*(1+cos(alpha));
    l = 1;
    
    % Study the flow in function of Q_dot
    
    % Define some usefull constant value 
    K       = sqrt(2*k/(k-1) * R/Mmol * (1 - (Pe/Pc)^((k-1)/k))); % V_e = k * sqrt(Tc)
    aa      = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v;
    cc      = cp_v*((T-Pe*Ae)/K)^2;
    Q_min   = sqrt(4*aa*cc);
    
    GAMMA = sqrt(k*(2/(k+1))^((k+1)/(k-1)));
    
    % f = @(x1, x2) interplinear(H2O.enthalpy, [x1 x2]);
    % F0 = @(Tc) Tc -(Pc*At/(Q_dot/((f(Tc, Pc*1e-6)-f(298.15, Pc*1e-6))*1e3))*GAMMA)^2*Mmol/R;
    % Tc = fzero(F0, 500)
    % m_dot = (T-Pe*Ae)/K/sqrt(Tc)
    % 
    % 
    % F0 = @(Tc) T- Q_dot/((f(Tc, Pc*1e-6)-f(298.15, Pc*1e-6))*1e3)*l*sqrt(Tc) - Pe*Ae;
    % Tc = fzero(F0,500)
    % m_dot = (T-Pe*Ae)/K/sqrt(Tc)
        
    % mass flow rate
    m_dot = (Q_dot + sqrt(Q_dot^2 - Q_min^2))/(2*aa)  % LOW TEMP BRANCH
    % m_dot = (Q_dot - sqrt(Q_dot^2 - Q_min^2))/(2*aa)  % HIGH TEMP BRANCH
    
    % Combustion Chamber Temperature
    Tc = ((T-Pe*Ae)/(m_dot *K))^2
    
    % Exit Velocity
    v_e = K *sqrt(Tc)
    
    % Exit Temperature
    T_e = Tc*(Pe/Pc)^((k-1)/k)
    
    % Specific Impulse
    Isp = T/(m_dot*g)
    
    % Injection Area
    dP_inj = dP_inj_perc*Pc
    Ainj = m_dot/(Cd*sqrt(2*dP_inj*rho_l));
    
    % dP_inj= 0.5*k_inj*m_dot^2/Ainj^2/rho_l;
    dP_pipe = 0.5*k_pipe*m_dot^2/Apipe^2/rho_l
    dP_v = 0.5*k_v*m_dot^2/Av^2/rho_l
    
    P_tank = Pc + dP_inj + dP_pipe + dP_v
    
    K_tot = 0.5/rho_l*(k_inj/Ainj^2 + k_pipe/Apipe^2 + k_v/Av^2);
    dP_tot = K_tot*m_dot^2;
    
    %% MONTECARLO
    
    % Definizione parametri
    tol = 1e-12;         % Tolleranza convergenza
        
    Apipe_std = (sqrt(Apipe/pi)+1e-6)^2*pi-Apipe;
    Apipe_avg = Apipe;
    Av_std = (sqrt(Av/pi)+1e-6)^2*pi-Av;
    Av_avg = Av;
    Ainj_std = (sqrt(Ainj/pi)+1e-6)^2*pi-Ainj;
    Ainj_avg = Ainj;
    At_std = (sqrt(At/pi)+1e-6)^2*pi-At;
    At_avg = At;
    Ae_std = (sqrt(Ae/pi)+1e-6)^2*pi-Ae;
    Ae_avg = Ae;
    Q_dot_std = Q_dot*1/100/3;
    Q_dot_avg = Q_dot;
    
    % Apipe_std = 0;
    % Av_std = 0;
    % Ainj_std = 0;
    % At_std = 0;
    % Ae_std = 0;
    % Q_dot_std = 0;
    
    Apipe_val = normrnd(Apipe_avg, Apipe_std, [N,1]);
    Av_val = normrnd(Av_avg, Av_std, [N,1]);
    Ainj_val = normrnd(Ainj_avg, Ainj_std, [N,1]);
    At_val = normrnd(At_avg, At_std, [N,1]);
    Ae_val = normrnd(Ae_avg, Ae_std, [N,1]);
    Q_dot_val =  normrnd(Q_dot_avg, Q_dot_std, [N, 1]);
    
    K_tot_val = 0.5/rho_l.*(k_inj./Ainj_val.^2 + k_pipe./Apipe_val.^2 + k_v./Av_val.^2);
    % dP_tot = K_tot*m_dot^2;
    
    % f = @(x1, x2) interplinear(H2O.enthalpy, [x1 x2]);
    % 
    % GAMMA = sqrt(k*(2/(k+1))^((k+1)/(k-1)));
    
    
    T_val = zeros(N,1);
    Pc_val = zeros(N,1);
    Pe_val = zeros(N,1);
    m_dot_val = zeros(N,1);
    ve_val = zeros(N,1);
    Tc_val = zeros(N,1);
    m_prop = zeros(N,1);
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
    m_prop_avg = zeros(N+1,1);
    m_prop_std = zeros(N+1,1);
    Isp_avg = zeros(N+1,1);
    Isp_std = zeros(N+1,1);
    
    
    unstart = 0;
    for i = 1:N
        Apipe = Apipe_val(i);
        Av    = Av_val(i);
        Ainj  = Ainj_val(i);
        At    = At_val(i);
        Ae    = Ae_val(i);
        K_tot = K_tot_val(i); % dP_tot = K_tot*m_dot^2;
        Q_dot = Q_dot_val(i);
    
        % it = 0;
        % err = 1;
        % m_dot =sqrt((P_tank-Pc)/K_tot);
        % 
        % while(err > tol && it < 1000)
        %     m_dot_old = m_dot;
        %     Tc_old = Tc;
        %     Pc_old = Pc;
        % 
        %     it = it + 1;
        %     Tc = (Pc*At/m_dot*GAMMA)^2*Mmol/R;
        %     m_dot = Q_dot/(f(Tc, Pc*1e-6)-f(298.15, Pc*1e-6))/1e3; 
        %     Pc = P_tank - K_tot*m_dot^2;
        %     if it > 999
        %         disp('It did not converge')
        %     end   
        %     err = abs(m_dot_old - m_dot)/m_dot+abs(Tc_old - Tc)/Tc+abs(Pc_old - Pc)/Pc;
        % end
        %
        % Pc = @(m_dot) P_tank - K_tot*m_dot.^2;
        % Tc = @(m_dot) ((P_tank - K_tot*m_dot.^2).*At./m_dot*GAMMA).^2*Mmol/R;
        % F0 = @(m_dot) m_dot*(f(Tc(m_dot), Pc(m_dot)*1e-6)-f(298.15, Pc(m_dot)*1e-6))*1e3 - Q_dot;
        % 
        % m_dot = fzero(F0, [1e-6 10e-6])
        % Tc = Tc(m_dot)
        % Pc = Pc(m_dot)
        %
        % m_dot = linspace(1e-6,5e-6,1000);
        % for i = 1:length(m_dot)
        %     try
        %     fff(i) = F0(m_dot(i));
        %     Pcc(i) = Pc(m_dot(i));
        %     Tcc(i) = Tc(m_dot(i));
        %     catch
        %         fff(i) = 0;
        %         Pcc(i) = 0;
        %         Tcc(i) = 0;
        %     end
        % end
        % subplot(1,3,1)
        % plot(m_dot,fff)
        % subplot(1,3,2)
        % plot(m_dot,Pcc)
        % subplot(1,3,3)
        % plot(m_dot,Tcc)
    
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
        p_ratio = xv;
        K_ratio2= k*(2/(k+1))^((k+1)/(k-1))*Mmol/R;
    
        AA = cp_v*K_ratio2*At^2*K_tot^2;
        BB = 0;
        CC = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v -2*P_tank*K_tot*cp_v*K_ratio2*At^2;
        DD = - Q_dot;
        EE = cp_v*K_ratio2*P_tank^2*At^2;
    
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
        %
    
        c1 = xmin;      % LOW TEMP BRANCH
        c2 = 1e-4;
        % c1 = 1e-7;  % HIGH TEMP BRANCH
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
        xv = x;
        % iter = 0;
        % while (iter< 15 && err> toll)
        %    dfx = dfun(xv);
        %    if dfx == 0
        %       error(' Arresto per azzeramento di dfun');
        %    else
        %       xn = xv - fun(xv)/dfx;
        %       err = abs(fun(xn));
        %       iter = iter+1;
        %       xv = xn;
        %    end
        % end
        m_dot = xv;
    
        Pc = P_tank- K_tot*m_dot^2;
        Tc = (Q_dot-m_dot*(cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v))/m_dot/cp_v;
        unstart = unstart + (Tc < T_eb);
        Tc = Tc*(Tc > T_eb);
        Pe = Pc*p_ratio;
        ve = sqrt(2*k/(k-1)*R/Mmol*Tc*(1-(Pe/Pc)^((k-1)/k)));
    
        T_val(i)=(m_dot*l*ve+Pe*Ae)*(Tc > T_eb);
        Pc_val(i)=Pc*(Tc > T_eb);
        Pe_val(i)=Pe*(Tc > T_eb);
        m_dot_val(i)=m_dot*(Tc > T_eb);
        ve_val(i)=ve*(Tc > T_eb);
        Isp_val(i) = T_val(i)/m_dot_val(i)/g;
        Tc_val(i) = Tc;
    
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
    
    unstart_design(iii) = unstart/N*100;
    T_avg_design(iii) = T_avg(end);
    T_std_design(iii) = T_std(end);
    Pc_avg_design(iii) = Pc_avg(end);
    Pc_std_design(iii) = Pc_std(end);
    Pe_avg_design(iii) = Pe_avg(end);
    Pe_std_design(iii) = Pe_std(end);
    m_dot_avg_design(iii) = m_dot_avg(end);
    m_dot_std_design(iii) = m_dot_std(end);
    ve_avg_design(iii) = ve_avg(end);
    ve_std_design(iii) = ve_std(end);
    Tc_avg_design(iii) = Tc_avg(end);
    Tc_std_design(iii) = Tc_std(end);
    Isp_avg_design(iii) = Isp_avg(end);
    Isp_std_design(iii) = Isp_std(end);
end

%% PLOTS
figure(1)
plot(Q_design,unstart_design)

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

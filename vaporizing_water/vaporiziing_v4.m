%% Vaporizing Liquid Microthrusters V4
%
% -------------------------------------------------------------------------
clear
clc
close all

%% DATA
eps     = 9;            % Rapporto tra le aree (ugello convergente-divergente)  [adim]
k       = 1.327;        % Rapporto dei calori specifici                         [adim]
g       = 9.81;         % Gravity acceleration                                  [m/s^2]
Pc      = 230000;       % Evaporation chamber Pressure                          [Pa]
R       = 8.314;        % Gas Constant                                          [J/mol*K]
Mmol    = 0.018;        % Water Molar Mass                                      [kg/mol]
rho_l   = 999;          % Density liquid water                                  [kg/m^3]
rho_v   = 0.6;          % Density gas water                                     [kg/m^3]
cp_l    = 4.1838e3;     % Specific Heat liquid water                            [J/kgK]
cp_v    = 2.0256e3;     % Specific Heat gas water                               [J/kgK]
lambda  = 2260e+3;      % Water latent evaporation heat                         [J/kg]
T_eb    = 124 + 273;    % Temperature of ebollition                             [K]
T_inj   = 25 + 273.15;  % temperature of Injection                              [K]
alpha   = deg2rad(30);  %                                                       [deg]
beta    = deg2rad(45);  %                                                       [deg]
Cd      = 0.7;          % Drag coefficient                                      [adim]
dP_inj  = 0.1;          % Pressure lost for injection   !!!CHECK!!!
T       = 5e-3;         % Wanted Thrust                                         [N]
dv      = 0.5;          % Wanted Delta v                                        [m/s]
mass    = 4;            % S/C mass                                              [kg]
N       = 100;          % Number of sample for Q_dot                            [adim]
k_gas   = 1.67;         % He 
R_gas = 2077.3;         %J/kgK

%% Calculate the Pressure Ratio (Pe/Pc)

% Equation
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

%% Compute throat and exit area 

CT = k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k)))+Pe/Pc*eps;
At = T/Pc/CT;
r_t = sqrt(At/pi);
Ae = eps*At;
r_e = sqrt(Ae/pi);

%% Dimensionate nozzle
eps_conv = 1/0.2 * sqrt((2/(k+1)*(1+(k-1)/2*0.2^2))^((k+1)/(k-1)));
Ac = eps_conv * At;
r_c = sqrt(Ac/pi);

L_conv = 1/2*(2*(r_c-r_t)/tan(beta)); % ipotesi che Ae è uagule a quwlla di ingresso nel convergente
L_div = 1/2*(2*(r_e-r_t)/tan(alpha));
L_nozzle = L_conv + L_div;
l = 1/2*(1+cos(alpha));         % 2D losses coefficientl

%% Pipes geometry from literature
D_h = 500e-6;
L = 1000e-6;             
mu = 1.137e-3;
A_pipe = pi*(D_h/2)^2;

%% Study the flow in function of Q_dot

% Define some usefull constant value 
K       = sqrt(2*k/(k-1) * R/Mmol * (1 - (Pe/Pc)^((k-1)/k))); % V_e = k * sqrt(T_c)
dd      = cp_l * (T_eb-T_inj) + lambda - T_eb*cp_v;
ee      = cp_v*((T-Pe*Ae)/(l*K))^2;
Q_min   = sqrt(4*dd*ee);
k_inj = 1/Cd^2; % Injection losses with Cd = 0.7

dt = dv * mass / T; % time of burn with max mass sat = 4 (max possible)

m_dot   = zeros(N,2);
T_c     = zeros(N,2);
v_e     = zeros(N,2);
T_e     = zeros(N,2);
Isp     = zeros(N,2);
A_inj   = zeros(N,2);
A_val   = zeros(N,2);
A_reg   = zeros(N,2);
v = zeros(N,2);
Re = zeros(N,2);
P_tank = zeros(N,2);
m_fuel = zeros(N,2);
V_fuel = zeros(N,2);
P_gas_i = zeros(N,2);
M_gas = zeros(N,2);
V_gas = zeros(N,2);
p_drop_injection = zeros(N,1);
P_plenum = zeros(N,1);
p_drop_pipes = zeros(N,2);
p_drop_valve = zeros(N,2);
p_drop_reg = zeros(N,1);
Q_dot = linspace(Q_min,13,N);


for i = 1:1
    Q_tot = 13;
    % mass flow rate
    m_dot(i,1) = (Q_dot(i) + sqrt(Q_dot(i)^2 - Q_min^2))/(2*dd);
    m_dot(i,2) = (Q_dot(i) - sqrt(Q_dot(i)^2 - Q_min^2))/(2*dd);
    
    % Combustion Chamber Temperature
    T_c(i,1) = ((T-Pe*Ae)/(m_dot(i,1) *K))^2;
    T_c(i,2) = ((T-Pe*Ae)/(m_dot(i,2) *K))^2;
    
    % Exit Velocity
    v_e(i,1) = K *sqrt(T_c(i,1));
    v_e(i,2) = K *sqrt(T_c(i,2));

    % Exit Temperature
    T_e(i,1) = T_c(i,1)*(Pe/Pc)^((k-1)/k);
    T_e(i,2) = T_c(i,2)*(Pe/Pc)^((k-1)/k);
    
    % Specific Impulse
    Isp(i,1) = T/(m_dot(i,1)*g);
    Isp(i,2) = T/(m_dot(i,2)*g);

    % Injection Area
    A_inj(i,1) = m_dot(i,1)/(Cd*sqrt(2*dP_inj*Pc*rho_l));
    A_inj(i,2) = m_dot(i,2)/(Cd*sqrt(2*dP_inj*Pc*rho_l));

    % Injection Losses 
    p_drop_injection(i) = 0.1*Pc; % ipotizzato 10% di perdita di pressione
    P_plenum(i) = Pc + p_drop_injection(i); % tale pressione è garantita dal regolatore   
   
    % Valve losses
    p_drop_valve(i) = 0.04*(Pc+p_drop_injection(i)); % da paper thesis 0.03, da paper hydrazine 0.04
    
    % Pipes losses
    v(i,1) = m_dot(i,1)/(A_pipe*rho_l);
    v(i,2) = m_dot(i,2)/(A_pipe*rho_l);
    Re(i,1) = rho_l.*v(i,1)*D_h/mu;
    Re(i,2) = rho_l.*v(i,2)*D_h/mu;
    p_drop_pipes(i,1) = 64./Re(i,1) * (L/D_h)*rho_l*mu^2/2;
    p_drop_pipes(i,2) = 64./Re(i,2) * (L/D_h)*rho_l*mu^2/2;
    
    % Tank Pressure
    P_tank(i, 1) = P_plenum(i) + p_drop_pipes(i,1) + p_drop_valve(i);
    P_tank(i, 2) = P_plenum(i) + p_drop_pipes(i,2) + p_drop_valve(i);

    % Mass propellant budget
    m_fuel(i,1) = m_dot(i,1)*dt*1.055;
    m_fuel(i,2) = m_dot(i,2)*dt*1.055; % 5.5% di margine dalle ECSS
    V_fuel(i,1) = m_fuel(i,1)/rho_l*1.1; % 1% di margine dalle ECSS
    V_fuel(i,2) = m_fuel(i,2)/rho_l*1.1;

    
    % Option with pressurant 
%     P_gas_i(i,1) = 5*P_tank(i, 1); % scelta
%     P_gas_i(i,2) = 5*P_tank(i, 2); 
%     M_gas(i,1) = P_tank(i,1)*V_fuel(i,1)/(R_gas*T_inj)*k_gas/(1+(P_tank(i,1)/P_gas_i(i,1)));
%     M_gas(i,2) = P_tank(i,2)*V_fuel(i,2)/(R_gas*T_inj)*k_gas/(1+(P_tank(i,2)/P_gas_i(i,2)));
%     V_gas(i,1) = M_gas(i,1)*R_gas*T_inj/P_gas_i;
%     V_gas(i,2) = M_gas(i,2)*R_gas*T_inj/P_gas_i;
    
end

%% Temperature, Pressure and Mach inside the nozzle
T_nozzle = [];
MACH =[];
P_nozzle = [];

% sampling points
points_nozzle = linspace(0,L_nozzle,2*N);
radius_nozzle_convergent = linspace(r_c,r_t,N);
radius_nozzle_divergent  = linspace (r_t,r_e,N);

% Calculate area at sampled points
Area_nozzle_convergent = pi * radius_nozzle_convergent.^2;
Area_nozzle_divergent = pi * radius_nozzle_divergent.^2;
Area_nozzle_vect = [Area_nozzle_convergent,Area_nozzle_divergent];
Area_ratio_vect = Area_nozzle_vect/At;  % Area Ratio

% Calculate Mach at different sampled points

options = optimset('TolX', 1e-8, 'Tolfun', 1e-8);

for i=[1:length(Area_ratio_vect)]
    A_ratio = Area_ratio_vect(i);
    fun_Mach = @(M) (1/M)*((2/(k+1))*(1 + ((k-1)/2)*M^2))^((k+1)/(2*(k-1))) - A_ratio;
    
    if i<=N
        mach = fsolve(fun_Mach,0.3);
        MACH = [MACH,mach];
    else
        mach = fsolve(fun_Mach,8);
        MACH = [MACH,mach];
    end
end

% Calculate temperature at sampled times
T_tot_vect = T_c * (1 + (k-1)/2 * MACH(1)^2);

T_tot       = T_tot_vect (50,2);
T_exit      = T_e (50,2);

for i=[1:length(Area_ratio_vect)]
    mach = MACH(i);
    denom = (1 + (k-1)/2*mach^2);
    temp = T_tot / denom;
    T_nozzle =[T_nozzle,temp];
end

% Calculate pressure at sampled times
P_nozzle = Pe .* (T_nozzle ./ T_exit) .^(k/(k-1));


%% PLOTS

figure
% Mass flow rate 
subplot(2,4,1)
plot(Q_dot,m_dot(:,1), "LineWidth",1)
hold on
plot(Q_dot,m_dot(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Mass flow Rate [kg/s]');
legend('Solution #1','Solution #2');

% Combustion Chamber Temperature
subplot(2,4,2)
plot(Q_dot,T_c(:,1), "LineWidth",1)
hold on
plot(Q_dot,T_c(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Combustion Chamber Temperature [K]');
legend('Solution #1','Solution #2');

% Exit Temperature
subplot(2,4,3)
plot(Q_dot,T_e(:,1), "LineWidth",1)
hold on
plot(Q_dot,T_e(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Exit Temperature [K]');
legend('Solution #1','Solution #2');

% Exit Velocity
subplot(2,4,4)
plot(Q_dot,v_e(:,1), "LineWidth",1)
hold on
plot(Q_dot,v_e(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Exit Velocity [m/s]');
legend('Solution #1','Solution #2');

% Specific Impulse
subplot(2,4,5)
plot(Q_dot,Isp(:,1), "LineWidth",1)
hold on
plot(Q_dot,Isp(:,2),"LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Specific Impulse [s]');
legend('Solution #1','Solution #2');

% Injection Area
subplot(2,4,6)
plot(Q_dot,A_inj(:,1), "LineWidth",1)
hold on
plot(Q_dot,A_inj(:,2), "LineWidth",1)
grid on;
xlabel('Thermal Power [W]');
ylabel('Injection Area [m^2]');
legend('Solution #1','Solution #2');


%% Plot Nozzle 2D and 3D

% Nozzle
% figure


Np = 100;       % Number of points for the plot
Ntheta = 50;    % Number of points for the rotation
% Convergent
x_conv = linspace(0, L_conv, Np);
r_conv = linspace(r_c, r_t, Np);

% Divergent
x_div = linspace(L_conv, L_conv + L_div, Np);
r_div = linspace(r_t, r_e, Np);

% Merge the sections
x_total = [x_conv, x_div];
r_total = [r_conv, r_div];

% Conical supposition
theta = linspace(0, 2*pi, Ntheta);  % Angle from 0 a 360°

% 3D mesh
[X, THETA] = meshgrid(x_total, theta);
R = repmat(r_total, length(theta), 1); % estende i raggi
Y = R .* cos(THETA);
Z = R .* sin(THETA);

% Plot the profile (fill + plot)
figure;
subplot(1,2,1)
hold on;
fill([x_total, fliplr(x_total)], [r_total, -fliplr(r_total)], [0.8, 0.8, 1]); % colore azzurrino
plot(x_total, r_total, 'b', 'LineWidth', 1.5); % profilo superiore
plot(x_total, -r_total, 'b', 'LineWidth', 1.5); % profilo inferiore
xlabel('Length [m]');
ylabel('Radiuns [m]');
title('Nozzle 2D Profile');
axis equal;
grid on;

% Plot 3D
subplot(1,2,2)
surf(X, Y, Z, 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none');
xlabel('Lunghezza [m]');
ylabel('Asse Y [m]');
zlabel('Asse Z [m]');
title('3D Nozzle');
axis equal;
grid on;
camlight; lighting gouraud; % Effetti di luce per renderlo più cute <3 <3 <3 <3

%% plot properties along the nozzle

% Mach number
[X, Y] = meshgrid(x_total, [-1, 1]); % 2 righe: parte superiore e inferiore
Y = Y .* r_total;                   % scala in base al raggio locale
C = repmat(MACH, 2, 1);      % matrice dei valori di Mach (2 righe uguali)

figure;
subplot(1,3,1)
hold on;
s = surf(X, Y, zeros(size(X)), C, 'EdgeColor', 'none'); % z=0 per renderlo 2D
view(2);                        % vista 2D
colormap(jet);                 % scegli la colormap (es. jet, parula, etc.)
colorbar;                      % aggiunge la barra del colore
xlabel('Length [m]');
ylabel('Radius [m]');
title('Nozzle 2D Profile Colored by Mach');
axis equal;
grid on;

% Temperature
[X, Y] = meshgrid(x_total, [-1, 1]); % 2 righe: parte superiore e inferiore
Y = Y .* r_total;                   % scala in base al raggio locale
C = repmat(T_nozzle, 2, 1);      % matrice dei valori di Mach (2 righe uguali)

subplot(1,3,2)
hold on;
s = surf(X, Y, zeros(size(X)), C, 'EdgeColor', 'none'); % z=0 per renderlo 2D
view(2);                        % vista 2D
colormap(jet);                 % scegli la colormap (es. jet, parula, etc.)
colorbar;                      % aggiunge la barra del colore
xlabel('Length [m]');
ylabel('Radius [m]');
title('Nozzle 2D Profile Colored by Temperature');
axis equal;
grid on;

% Pressure
[X, Y] = meshgrid(x_total, [-1, 1]); % 2 righe: parte superiore e inferiore
Y = Y .* r_total;                   % scala in base al raggio locale
C = repmat(P_nozzle, 2, 1);      % matrice dei valori di Mach (2 righe uguali)

subplot(1,3,3)
hold on;
s = surf(X, Y, zeros(size(X)), C, 'EdgeColor', 'none'); % z=0 per renderlo 2D
view(2);                        % vista 2D
colormap(jet);                 % scegli la colormap (es. jet, parula, etc.)
colorbar;                      % aggiunge la barra del colore
xlabel('Length [m]');
ylabel('Radius [m]');
title('Nozzle 2D Profile Colored by Pressure');
axis equal;
grid on;
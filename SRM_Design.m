clear
clc
close all

codes = 2438:2446;
trigger_action = 0.05; % 5%
f = 1000; % hz
n = [1,9999];
web = 30; % mm
L = 290; % mm
D = 160; % mm
frac_ap = 0.68;
frac_al = 0.18;
frac_htpb = 0.14;
rho_ap = 1.95; % g/cm3
rho_al = 2.7; % g/cm3
rho_htpb = 0.92; % g/cm3
Dt = [28.80, 25.26, 21.81]; % mm

pbarStruct = load('tracesbar1.mat'); % bar

p_low = 0;
p_mid = 0;
p_high = 0;

for i = n(1):min(length(codes), n(2))
	name = join(["pbarStruct.pbar",string(codes(i))],"");
	p_array = eval(name);
	p0 = p_array(1,:);
	%p_array = p_array - p0;
	if i ~= 7
	p_low(1:length(p_array(:,1)),i) = p_array(:,1);
	p_mid(1:length(p_array(:,2)),i) = p_array(:,2);
	p_high(1:length(p_array(:,3)),i) = p_array(:,3);
	else
	p_low(1:length(p_array(:,1)),i) = p_array(:,1);
	p_mid(1:length(p_array(:,2)),i) = p_array(:,3);
	p_high(1:length(p_array(:,3)),i) = p_array(:,2);
	end

	p_max(i,1) = max(p_low(:,i));
	p_max(i,2) = max(p_mid(:,i));
	p_max(i,3) = max(p_high(:,i));

	% p_action = trigger_action * p_max;
	p_action = trigger_action * (p_max - p0) + p0;

	t_action(i,1,:) = [find(p_low(1:floor(end/2),i) <= p_action(i,1), 1, "last"), find(p_low(:,i) >= p_action(i,1), 1, "last")];
	t_action(i,2,:) = [find(p_mid(1:floor(end/2),i) <= p_action(i,2), 1, "last"), find(p_mid(:,i) >= p_action(i,2), 1, "last")];
	t_action(i,3,:) = [find(p_high(1:floor(end/2),i) <= p_action(i,3), 1, "last"), find(p_high(:,i) >= p_action(i,3), 1, "last")];

	p_ref(i,1) = sum(p_low(t_action(i,1,1):t_action(i,1,2),i)) / 2 / (t_action(i,1,2) - t_action(i,1,1));
	p_ref(i,2) = sum(p_mid(t_action(i,2,1):t_action(i,2,2),i)) / 2 / (t_action(i,2,2) - t_action(i,2,1));
	p_ref(i,3) = sum(p_high(t_action(i,3,1):t_action(i,3,2),i)) / 2 / (t_action(i,3,2) - t_action(i,3,1));

	t_burn(i,1,:) = [find(p_low(1:floor(end/2),i) <= p_ref(i,1), 1, "last"), find(p_low(:,i) >= p_ref(i,1), 1, "last")];
	t_burn(i,2,:) = [find(p_mid(1:floor(end/2),i) <= p_ref(i,2), 1, "last"), find(p_mid(:,i) >= p_ref(i,2), 1, "last")];
	t_burn(i,3,:) = [find(p_high(1:floor(end/2),i) <= p_ref(i,3), 1, "last"), find(p_high(:,i) >= p_ref(i,3), 1, "last")];

	p_eff(i,1) = sum(p_low(t_burn(i,1,1):t_burn(i,1,2),i)) / (t_burn(i,1,2) - t_burn(i,1,1));
	p_eff(i,2) = sum(p_mid(t_burn(i,2,1):t_burn(i,2,2),i)) / (t_burn(i,2,2) - t_burn(i,2,1));
	p_eff(i,3) = sum(p_high(t_burn(i,3,1):t_burn(i,3,2),i)) / (t_burn(i,3,2) - t_burn(i,3,1));
end

rb = web ./ (t_burn(:,:,2) - t_burn(:,:,1)) * f;

rb_calc = [rb(:,1); rb(:,2); rb(:,3)];
p_calc = [p_eff(:,1); p_eff(:,2); p_eff(:,3)];

[a_avg, a_std, n_avg, n_std, R2] = Uncertainty(p_calc, rb_calc);

rho = 1/ (frac_ap/rho_ap + frac_al/rho_al + frac_htpb/rho_htpb) * 1e3; % kg/m3

clearvars -except a_avg a_std n_avg n_std rho

%% NOMINAL OPERATION

% Dati del problema
eps = 2.6      % Rapporto tra le aree (ugello convergente-divergente)
k = 1.21;       % Rapporto dei calori specifici
Pc = 5e+5    % Pressione in camera di combustione (Pa) DA GUARDARE DA LETTERATURA
lll=0.912%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DA VEDERE FATTO CON CHAT, RE_d=723
%FATTO CONSIDERANDO B=inf che è diversissimo ma sperimentalmente si
%verifica che cambia poco
alpha = 30;
lambda = (1 + cosd(alpha))/2;
% lambda=1;

T = 5e-3; 
DeltaV = 0.5; % [m/s]
Mass = 4; % [kg]

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
Pe = xn * Pc

% Parametri aggiuntivi (devono essere definiti prima dell'uso)
R = 8.314;        % Costante specifica del gas (J/mol*K), da definire correttamente per il propellente
Mmol = 0.028;   % Massa molare del gas (kg/mol), da definire correttamente
Tc = 3070;      % Temperatura in camera di combustione (K), da definire

% Calcolo della velocità di scarico
ve = sqrt(2*k/(k-1) * R/Mmol * Tc * (1 - (Pe/Pc)^((k-1)/k) ) )
CT = k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-(Pe/Pc)^((k-1)/k)))*lll+Pe/Pc*eps;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
At = T/Pc/CT;
Ae = eps*At;
m_dot = (T-Pe*Ae)/(ve*lambda)/lll%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_avg = a_avg*10^(-3-5*n_avg);
a_std = a_std*10^(-3-5*n_avg);
Ab = m_dot/(Pc^n_avg*a_avg*rho);
h = a_avg*Pc^n_avg*DeltaV/T*Mass;

T_g = Tc/(1+ (k-1)/2)

%% MONTECARLO
N = 20000;

% manca il valore esatto di incertezza
At_std = (sqrt(At/pi)+0.001*1e-3)^2*pi-At;
At_avg = At;
Ae_std = (sqrt(Ae/pi)+0.001*1e-3)^2*pi-Ae;
Ae_avg = Ae;

% At_std = 0;
% Ae_std = 0;
% a_std = 0;
% n_std = 0;

At_val = normrnd(At_avg, At_std, [N,1]);
Ae_val = normrnd(Ae_avg, Ae_std, [N,1]);
a_val =  normrnd(a_avg, a_std, [N, 1]);
n_val =  normrnd(n_avg, n_std, [N, 1]);

T_val=zeros(N,1);
Pc_val=zeros(N,1);
Pe_val=zeros(N,1);
m_dot_val=zeros(N,1);
ve_val=zeros(N,1);
h_val=zeros(N,1);
m_prop =zeros(N,1);

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
h_avg = zeros(N+1,1);
h_std = zeros(N+1,1);
m_prop_avg = zeros(N+1,1);
m_prop_std = zeros(N+1,1);

for i=1:N
    At_mc=At_val(i);
    Ae_mc=Ae_val(i);
    a_mc=a_val(i);
    n_mc=n_val(i);
    eps_mc = Ae_mc/At_mc;

    % PRECACLOLARE LE COSTANTI DI K
    % eq1 = T/(Pc*At)==k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1)))...
    %     *sqrt(1-(x)^((k-1)/k))+eps*x;
    % eq2 = T == m_dot*ve + Pe*Ae;
    % eq3 = ve == sqrt(2*k/(k-1) * R/Mmol * Tc * (1 - (x)^((k-1)/k)));
    % eq4 = m_dot == rho*Ab*Pc^n*a;
    % eq5 = 0 == -1/eps + ((k+1)/2)^(1/(k-1)) * x^(1/k) * sqrt((k+1)/(k-1) * (1-(x)^((k-1)/k)));
	% eq1 = rho*Ab/At*Pc^(n-1)*a*sqrt(R/Mmol*Tc) == sqrt(k*(2/(k+1))^((k+1)/(k-1)));
	% eq2 = 1/eps == ((k+1)/2)^(1/(k-1))*x^(1/k)*sqrt((k+1)/(k-1)*(1-x^((k-1)/k)));
    % sol=vpasolve([eq1, eq2], [Pc,Pe]);
	% Pc = double(sol.Pc);
	% Pe = double(sol.Pe);
    % x = Pe/Pc;

    K = sqrt(k*(2/(k+1))^((k+1)/(k-1))/R/Tc*Mmol)/rho/Ab;
    Pc = (K*At_mc/a_mc/lambda)^(1/(n_mc-1));
    
    aa = 1/eps_mc^2*((k+1)/2)^(2/(1-k))*(k-1)/(k+1);
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
    Pe = xv*Pc;

	m_dot = rho*Ab*Pc^n_mc*a_mc*lll;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ve = sqrt(2*k/(k-1) * R/Mmol * Tc * (1 - (Pe/Pc)^((k-1)/k)));

    T_val(i)=m_dot*lambda*ve+Pe*Ae_mc;
    Pc_val(i)=Pc;
    Pe_val(i)=Pe;
    m_dot_val(i)=m_dot;
    ve_val(i)=ve;

    h_val(i)= a_mc*Pc_val(i)^n_mc*DeltaV/T_val(i)*Mass;
    m_prop(i) = h_val(i)*Ab*rho;
    
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
    h_avg(i+1)=(h_avg(i)*(i-1)+h_val(i))/i;
    h_std(i+1)=std(h_val(1:i));
    m_prop_avg(i+1)=(m_prop_avg(i)*(i-1)+m_prop(i))/i;
    m_prop_std(i+1)=std(m_prop(1:i));

end

%%

subplot(2,6,1)
plot(T_avg)
subplot(2,6,7)
plot(T_std)
subplot(2,6,2)
plot(Pc_avg)
subplot(2,6,8)
plot(Pc_std)
subplot(2,6,3)
plot(Pe_avg)
subplot(2,6,9)
plot(Pe_std)
subplot(2,6,4)
plot(m_dot_avg)
subplot(2,6,10)
plot(m_dot_std)
subplot(2,6,5)
plot(ve_avg)
subplot(2,6,11)
plot(ve_std)
subplot(2,6,6)
plot(h_avg)
subplot(2,6,12)
plot(h_std)

%% Array Sizing

wall_th = 1e-3;
d_b = sqrt(Ab/pi)*2
d_t = sqrt(At/pi)*2
d_e = sqrt(Ae/pi)*2
d_max = max([d_b, d_e]);

N_side = round(sqrt((h+5*h_std(end))/1e-3))
N_thrusters = N_side^2

L_plate = (d_max + wall_th)*N_side + wall_th

r_b = d_b/2;
r_t = d_t/2;
r_e = d_e/2;
L_conv = (r_b-r_t)/tan(deg2rad(45))
L_div = (r_e-r_t)/tan(deg2rad(30))
L = L_conv+L_div

h_prop = (h+5*h_std(end))/N_thrusters
clear
clc

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
figure(1)
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

[a_avg, a_std, n_avg, n_std, R2] = Uncertainty(p_calc, rb_calc)

rho = 1/ (frac_ap/rho_ap + frac_al/rho_al + frac_htpb/rho_htpb) * 1e3; % kg/m3

clearvars -except a_avg a_std n_avg n_std rho

%%
R = 8.314;        % Costante specifica del gas (J/kg*K), da definire correttamente per il propellente
Mmol = 0.028;   % Massa molare del gas (kg/mol), da definire correttamente
Tc = 3000;      % Temperatura in camera di combustione (K), da definire
k = 1.4;
Ae =  2.673591011270363e-07;
Ab =  1.261940759218272e-06;

At_avg =  2.673591011270363e-08;
a_stdt = At_avg*1e-5;

N = 50;
alpha = 30;
lambda = (1 + cos(alpha))/2;

At_val = normrnd(At_avg, a_stdt, [N,1]);
a_val = normrnd(a_avg, a_std, [N, 1]);
n_val = normrnd(n_avg, n_std, [N, 1]);
a_val = a_val.*10.^(-3-5*n_val);

T_val=zeros(N,1);
Pc_val=zeros(N,1);
Pe_val=zeros(N,1);
m_dot_val=zeros(N,1);
ve_val=zeros(N,1);

for i=1:N
    
    syms T Pc Pe m_dot ve positive
    
    assumeAlso(T,"real");
    assumeAlso(Pc,"real");
    assumeAlso(Pe,"real");
    assumeAlso(m_dot,"real");
    assumeAlso(ve,"real");
    
    At=At_val(i);
    a=a_val(i);
    n=n_val(i);
    x = Pe/Pc;
    eps = Ae/At;
    % eq1 = T/(Pc*At)==k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1)))...
    %     *sqrt(1-(x)^((k-1)/k))+eps*x;
    % eq2 = T == m_dot*ve + Pe*Ae;
    % eq3 = ve == sqrt(2*k/(k-1) * R/Mmol * Tc * (1 - (x)^((k-1)/k)));
    % eq4 = m_dot == rho*Ab*Pc^n*a;
    % eq5 = 0 == -1/eps + ((k+1)/2)^(1/(k-1)) * x^(1/k) * sqrt((k+1)/(k-1) * (1-(x)^((k-1)/k)));
    
	eq1 = rho*Ab/At*Pc^(n-1)*a*sqrt(2*k/(k-1)*R/Mmol*Tc*(1-x^((k-1)/k))) == k*sqrt(2/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1-x^((k-1)/k)));
	eq2 = 1/eps == ((k+1)/2)^(1/(k-1))*x^(1/k)*sqrt((k+1)/(k-1)*(1-x^((k-1)/k)));
    sol=vpasolve([eq1, eq2], [Pc,Pe]);
	Pc = double(sol.Pc);
	Pe = double(sol.Pe);
    x = Pe/Pc;

	eq3 = m_dot == rho*Ab*Pc^n*a;
    m_dot = double(vpasolve(eq3, m_dot));

    eq4 = ve == sqrt(2*k/(k-1) * R/Mmol * Tc * (1 - (x)^((k-1)/k)));
    ve = double(vpasolve(eq4, ve));

    T_val(i)=m_dot*lambda*ve+Pe*Ae;
    Pc_val(i)=sol.Pc;
    Pe_val(i)=sol.Pe;
    m_dot_val(i)=m_dot;
    ve_val(i)=ve;

end
%%
T_val=abs(T_val);
Pc_val=abs(Pc_val);
Pe_val=abs(Pe_val);
m_dot_val=abs(m_dot_val);
ve_val=abs(ve_val);

for i=1:N
    T_avg(i)=mean(T_val(1:i));
    T_std(i)=std(T_val(1:i));
    Pc_avg(i)=mean(Pc_val(1:i));
    Pc_std(i)=std(Pc_val(1:i));
    Pe_avg(i)=mean(Pe_val(1:i));
    Pe_std(i)=std(Pe_val(1:i));
    m_dot_avg(i)=mean(m_dot_val(1:i));
    m_dot_std(i)=std(m_dot_val(1:i));
    ve_avg(i)=mean(ve_val(1:i));
    ve_std(i)=std(ve_val(1:i));
end

subplot(2,5,1)
plot(T_avg)
subplot(2,5,6)
plot(T_std)
subplot(2,5,2)
plot(Pc_avg)
subplot(2,5,7)
plot(Pc_std)
subplot(2,5,3)
plot(Pe_avg)
subplot(2,5,8)
plot(Pe_std)
subplot(2,5,4)
plot(m_dot_avg)
subplot(2,5,9)
plot(m_dot_std)
subplot(2,5,5)
plot(ve_avg)
subplot(2,5,10)
plot(ve_std)
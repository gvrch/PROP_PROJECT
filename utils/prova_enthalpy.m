close all
load('CH4_20_30MPa_300_625K_Test')
idx2 = find(CH4_20_30MPa_300_625K_Test.enthalpy.d{2}<1000);
std = scraping4real(298.15, 27, 27, 0, 'CH4', 'Enthalpy', 1);
idx = find(CH4_20_30MPa_300_625K_Test.enthalpy.d{1}==27);
figure
plot(CH4_20_30MPa_300_625K_Test.enthalpy.d{2}, CH4_20_30MPa_300_625K_Test.enthalpy.data(idx, :)-std, 'k');
p1 = polyfit(CH4_20_30MPa_300_625K_Test.enthalpy.d{2}(idx2), CH4_20_30MPa_300_625K_Test.enthalpy.data(idx, idx2), 1);
p3 = polyfit(CH4_20_30MPa_300_625K_Test.enthalpy.d{2}, CH4_20_30MPa_300_625K_Test.enthalpy.data(idx, :), 3);
x = 100:1:1000;
y1 = polyval(p1, x)-std;
y3 = polyval(p3, x)-std;
hold on
% plot(x, y1, 'b');
% plot(x, y3, 'g');
run("HGdata_CH4.m")

std1 = scraping4real(298.15, 0.101325, 0.101325, 0, 'CH4', 'Enthalpy', 1);
yform = (Enthalpy(CH4.a, CH4.b1, x, CH4.H0)-Enthalpy(CH4.a, CH4.b1, 298.15, CH4.H0))/(16.0425e-3);
plot(x, yform);

control = scraping4real(600, 0.101325, 0.101325, 0, 'CH4', 'Enthalpy', 1)-std1;
scatter(600, control, '*');


model = fittype('(f*x^4 + a*x^3 + b*x^2 + c*x + d + e*(1/x))*x');
fit_mod = fit(CH4_20_30MPa_300_625K_Test.enthalpy.d{2}(idx2)', CH4_20_30MPa_300_625K_Test.enthalpy.data(idx, idx2)'-std, model);
y_mod =@(x) x.*(fit_mod.f.*x.^4+ fit_mod.a.*x.^3 + fit_mod.b.*x.^2 + fit_mod.c.*x + fit_mod.d + fit_mod.e.*(1./x));
% plot(x, y_mod(x), 'm')


p4 = polyfit(CH4_20_30MPa_300_625K_Test.enthalpy.d{2}(idx2), CH4_20_30MPa_300_625K_Test.enthalpy.data(idx, idx2), 2);
y4 = polyval(p4, x)-std;
% plot(x, y4, 'g');
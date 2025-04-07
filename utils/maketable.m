
% 
% Pmin = 15;
% Pmax = 35; 
% Pinc = 0.1;
% Tspan = 400:1:1100;
% 
% C = [];
% Cd2=[];
% for ii = 1:length(Tspan)
% x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2', 'Entropy');
% C.entropy.data(:, ii)=x.data;
% Cd2=[Cd2, x.d{2}];
% end
% C.entropy.d{2} = Cd2;
% C.entropy.d{1} = x.d{1};


Pmin = 0.1;
Pmax = 0.4; 
Pinc = 0.01;
Tspan = 293:4:500;

W=[];
C = [];
Cd2=[];
for ii = 1:length(Tspan)
x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2O', 'Enthalpy');
C.enthalpy.data(:, ii)=x.data;
Cd2=[Cd2, x.d{2}];
end
C.enthalpy.d{2} = Cd2;
C.enthalpy.d{1} = x.d{1};



% for ii = 1:length(Tspan)
% x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2O', 'Enthalpy');
% W.enthalpy.data(:, ii)=x.data;
% W=[W, x.d{2}];
% end
% % C.enthalpy.d{2}=Cd2;
% % C.enthalpy.d{1} = x.d{1};
% % std = scraping4real(298.15, 0.101325, 0.101325, 0,'CH4', 'Enthalpy');
% % C.enthalpy.std = std.data;
% 
% H2O = W


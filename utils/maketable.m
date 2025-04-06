
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


Pmin = 0.05;
Pmax = 1;
Pinc = 0.01;
Pspan = Pmin:Pinc:Pmax;

Tmin = 275;
Tmax = 1450; 
Tinc = 5;
Tspan = Tmin:Tinc:Tmax;

W=[];
C = [];
Cd2=[];
for ii = 1:length(Pspan)
x=scraping4real(Pspan(ii), Tspan(1), Tspan(end), Tinc,'H2O', 'Enthalpy',0,1);
C.enthalpy.data(:, ii)=x.data;
Cd2=[Cd2, x.d{2}];
end
C.enthalpy.d{2} = Cd2;
offset = 0;
for i=1:length(x.d{1})
    if sum(x.d{1}(i) == Tspan)
        C.enthalpy.d{1}(i-offset) = x.d{1}(i);
    else
        offset = 1;
    end
end
%%
H2O = C;
save("H2O_data.mat", "H2O")



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


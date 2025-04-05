Pmin = 20;
Pmax = 40; 
Pinc = 0.1;
Tspan = 50:2:1000;
Cd2=[];
for ii = 1:length(Tspan)
x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2', 'Density');
C.data(:, ii)=x.data;
Cd2=[Cd2, x.d{2}];
end
C.d{2}=Cd2;
C.d{1} = x.d{1};
rhoH2tab = C;
clear C
Cd2=[];
for ii = 1:length(Tspan)
x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2', 'Visc');
C.data(:, ii)=x.data;
Cd2=[Cd2, x.d{2}];
end
C.d{2}=Cd2;
C.d{1} = x.d{1};
viscH2tab = C;
clear C
Cd2=[];
for ii = 1:length(Tspan)
x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2', 'Cond');
C.data(:, ii)=x.data;
Cd2=[Cd2, x.d{2}];
end
C.d{2}=Cd2;
C.d{1} = x.d{1};
condH2tab = C;
clear C
Cd2=[];
for ii = 1:length(Tspan)
x=scraping4real(Tspan(ii), Pmin, Pmax, Pinc,'H2', 'Cp');
C.data(:, ii)=x.data;
Cd2=[Cd2, x.d{2}];
end
C.d{2}=Cd2;
C.d{1} = x.d{1};
CpH2tab = C;
url ='https://webbook.nist.gov/cgi/fluid.cgi?T=99&PLow=0&PHigh=10&PInc=0.5&Digits=5&ID=C7782447&Action=Load&Type=IsoTherm&TUnit=K&PUnit=MPa&DUnit=mol%2Fl&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF';
A = webread(url);
b = strfind (A,'<tr><th>Temperature (K)</th><th>Pressure (MPa)</th><th>Density (mol/l)</th><th>Volume (l/mol)</th><th>Internal Energy (kJ/kg)</th><th>Enthalpy (kJ/kg)</th><th>Entropy (J/g*K)</th><th>Cv (J/g*K)</th><th>Cp (J/g*K)</th><th>Sound Spd. (m/s)</th><th>Joule-Thomson (K/MPa)</th><th>Viscosity (uPa*s)</th><th>Therm. Cond. (W/m*K)</th><th>Phase</th>');
e = strfind (A,'</table');
e = e(e>=b);
e = e(1);

A = extractBetween(A, b, e);
A = split(A, '<tr>');
A = A(5:end);
A = split(A, '<td align="right">');
A = erase(A, '</td>');
A = erase(A, '</tr>↵');
A = A(:, 2:end-1);
B.data = A(3:end);
B.P = cell2mat(A(:, 2));
B.T = cell2mat(A(1,1));
B.fields = {'Density', 'Volume','Internal Energy', 'Enthalpy', 'Entropy','Cv (J/g*K)','Cp (J/g*K)','c (m/s)', 'Joule-Thomson (K/MPa)','Viscosity (uPa*s)','Therm. Cond.' };


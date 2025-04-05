function [B, Hf] = scraping4real(T, Pmin, Pmax, Pinc, substance, property, nofields, isobar)
% id methane = C74828
% id oxy = C7782447
switch substance
    case 'O2'
        id = 'C7782447';
        Hf = 0;
    case 'CH4'
        id = 'C74828';
        Hf = -74.87/16.0425e-3;
    case 'H2'
        id = 'C1333740';
        Hf = 0;
    case 'H2O'
        id = 'C7732185'; %max 2000K
        Hf =-241.826/18.0153e-3;
    case 'CO'
        id = 'C630080'; %max 500K
        Hf = -110.53/28.0101e-3;
    case 'CO2'
        id = 'C124389'; %max 2000K
        Hf = -393.51/44.0095e-3;
    case 'C2H6'
        id = 'C74840'; %max 675K
        Hf = -84/30.069e-3;
    case 'C2H4'
        id = 'C74851'; %max 450K
        Hf = 52.47/28.0532e-3;
    case 'C3H8'
        id = 'C74986'; %max 650K
        Hf = -104.7/44.0956e-3;
    
    otherwise 
        error('substance id not present')
end
if(nargin<=7||~isobar)
url =strcat('https://webbook.nist.gov/cgi/fluid.cgi?','T=',num2str(T),'&PLow=',...
    num2str(Pmin), '&PHigh=', num2str(Pmax), '&PInc=', num2str(Pinc), ...
    '&Digits=5&ID=',id,'&Action=Load&Type=IsoTherm&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF');
else
    url =strcat('https://webbook.nist.gov/cgi/fluid.cgi?','P=',num2str(T),'&TLow=',...
    num2str(Pmin), '&THigh=', num2str(Pmax), '&TInc=', num2str(Pinc), ...
    '&Digits=5&ID=',id,'&Action=Load&Type=IsoBar&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF');

end
options = weboptions('Timeout', 600);
A = webread(url, options);
b = strfind (A,'<tr><th>Temperature (K)</th>');
e = strfind (A,'</table');
e = e(e>=b);
e = e(1);

A = extractBetween(A, b, e);
A = split(A, '<tr>');

A= A(contains(A, '<td'));
A = split(A, '<td align="right">');
A = erase(A, '</td>');
A = erase(A, '</tr>↵');
if(iscolumn(A))
    A=A';
end
% A = A(contains(A(:, end), state), :);

A = A(:, 2:end-1); %!!!
if(nargin<=7||~isobar)
    col = A(:, 2);
    val = A(1, 1);
else
    col = A(:, 1);
    val = A(1, 2);
end
% A=extractBetween(A, 1, p);

for ii = 1:size(A, 2)
p = max(strlength(A(:,ii)));
for jj = 1:length(A(:, ii))

    for kk = 1:p-strlength(A(jj,ii))
        A(jj,ii)=strcat(A(jj,ii), '0');
    end

end
end
[col, idx] = unique(col);
[idx, id] = sort(idx);
col = col(id);
A = A(idx, :);
if(~isempty(Pinc))
cond = (mod(str2num(cell2mat(col)),Pinc)==0) | (Pinc ==0);
A = A(cond,:);
end
B.data = A(:,3:end);

B.d{1} = str2num(cell2mat(col));
B.d{2} = str2num(cell2mat(val));
fields = {'Density', 'Volume','Internal Energy', 'Enthalpy', 'Entropy','Cv (J/g*K)','Cp (J/g*K)','c (m/s)', 'Joule-Thomson (K/MPa)','Viscosity (uPa*s)','Therm. Cond.' };
try
idx = contains(fields, property);
B.data = str2num(cell2mat(B.data(:, idx)));
catch
    try
        idx = contains(fields, property(2:end));
        B.data = str2num(cell2mat(B.data(:, idx)));
    catch
         error('check property syntax')
    end
end
if(nargin>6 && nofields)
    B = B.data;
end
end
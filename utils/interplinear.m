function [datainterp]=interplinear(table, value)
% value = [T, P]
datainterp = table.data;
for i = 1 : length(value)
    l = find(value(i)>=table.d{i});
    u = find(value(i)<table.d{i});
    if(isempty(u))
        u = length(table.d{i});
        l = length(table.d{i})-1;
        warning('INTERPLINEAR: extrapolating')
    end
    if(isempty(l))
        u = 2;
        l = 1;
        warning('INTERPLINEAR: extrapolating')
    end
    u = u(1);
    l = l(end);
    if(isrow(datainterp))
    datainterp=datainterp';
    end

    datainterp = datainterp(l, :)+(value(i)-table.d{i}(l))/(table.d{i}(u)-table.d{i}(l)) .* (datainterp(u, :)-datainterp(l, :)); 
    datainterp=datainterp';
end
if(isfield(table, 'std'))
datainterp = datainterp - table.std;
end
end

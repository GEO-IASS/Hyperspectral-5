function o = multilin_options(select,linsolve,indepp,freep,acc)
o.select = select; %0 = multilineair, 1 = lineair, 2 = AAM
o.linsolve = linsolve;
o.indepp = indepp;
o.freeP = freep;
o.accuracy = acc;

names1 = {'multilineair','inverslineair','semilineair','lineair','AAM ML','AAM','semibilineair','bilineair'};
markers1 = {'h','','^','s','v','d'};
names2 = {'real','free','strict','strong'};
markers2 = {'b','m','g','r'};
name1 = names1{2*select+linsolve+1};
name2 = names2{freep+2*indepp+1};
marker1 = markers1{2*select+linsolve+1};
marker2 = markers2{freep+2*indepp+1};
if select && linsolve
    name2 = '';
    marker2 = 'k';
end

o.name = [name2 ' ' name1 ' (' num2str(acc) ')'];
o.marker = [marker1 marker2];

o.dof = 3 + (~linsolve)*(1+3*~indepp);

end
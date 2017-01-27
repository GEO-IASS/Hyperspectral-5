function o = multilin_options(select,linsolve,indepp,freep,acc)
o.select = select; %0 = multilineair, 1 = lineair, 2 = AAM, 3 = bilineair
o.linsolve = linsolve;
o.indepp = indepp;
o.freeP = freep;
o.accuracy = acc;
names1 = {'multilineair','inverslineair','semilineair','lineair','AAM ML','AAM','semibilineair','bilineair'};
names2 = {'real','free','strict','strong'};
name1 = names1{2*select+linsolve+1};
name2 = names2{freep+2*indepp+1};
if select && linsolve
    name2 = '';
end

o.name = [name2 ' ' name1 ' (' num2str(acc) ')'];

o.dof = 4 + (~linsolve)*(1+3*~indepp);

end
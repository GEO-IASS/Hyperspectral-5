function o = multilin_options(linapprox,linsolve,indepp,freep,acc)
o.linapprox = linapprox;
o.linsolve = linsolve;
o.indepp = indepp;
o.freeP = freep;
o.accuracy = acc;
names1 = {'multilineair','semilineair','inverslineair','lineair'};
names2 = {'real','free','strict','strong'};
name1 = names1{linapprox+2*linsolve+1};
name2 = names2{freep+2*indepp+1};
if linapprox && linsolve
    name2 = '';
end

o.name = [name2 ' ' name1 ' (' num2str(acc) ')'];

o.dof = 4 + (~linsolve)*(1+3*~indepp);

end
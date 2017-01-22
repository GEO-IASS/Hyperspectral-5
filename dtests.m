function s = dtests
s.pmc = @PMC_default;
s.pc = @pixelcompare;
s.l = @linfix;
s.mc = @mean_correctness;
s.tabletest = @Table_test_Alina;
s.tablehandle = @table_show_result;
end

function PMC_default(points,hits)
    PMC(points,hits,1,points,0.1,['export_PMC_' num2str(points) '_' num2str(hits) '_' ])
end


function PMC(points,hits,maxP,zoom_points,zoomP,filename)
    file_for = filename; 
    close all;
    
    smallP = linspace(0,zoomP,zoom_points+1);
    smallP = smallP(1:end-1);
    P = [smallP linspace(zoomP,maxP,points)];
    [mc,sc,me,se,mp,sp,ma,sa,ml,sl] = mean_correctness(P,hits);
    
    figure;
    errorbar(P,mc,sc);
    xlabel('P');
    ylabel('correct geselcteerde materialen');
    axis([0 Inf 0 Inf])
    axes('position',[.60 .60 .25 .25])
    errorbar(smallP,mc(1:zoom_points),sc(1:zoom_points));
    axis([0 Inf 0 Inf])
    print([file_for 'corr' '.pdf'],'-dpdf')
   
    figure;
    errorbar(P,me,se);
    xlabel('P');
    ylabel('reconstructieerror');
    axis([0 Inf 0 Inf])
    axes('position',[.2 .60 .25 .25])
    errorbar(smallP,me(1:zoom_points),se(1:zoom_points));
    axis([0 Inf 0 Inf])
    print([file_for 'E' '.pdf'],'-dpdf')
    
    figure;
    errorbar(P,me,se);
    hold on;
    errorbar(P,ml,sl,'r')
    hold off;
    xlabel('P');
    ylabel('reconstructieerror');
    axis([0 Inf 0 Inf])
    axes('position',[.2 .60 .25 .25])
    errorbar(smallP,me(1:zoom_points),se(1:zoom_points));
    hold on;
    errorbar(smallP,ml(1:zoom_points),sl(1:zoom_points));
    hold off;
    axis([0 Inf 0 Inf])
    print([file_for 'El' '.pdf'],'-dpdf')
    
    figure;
    h1 = errorbar(P,me,se);
    hold on;
    h2 = errorbar(P,ml,sl,'r');
    hold off;
    set(get(h1,'Parent'), 'YScale', 'log')
    set(get(h2,'Parent'), 'YScale', 'log')
    axis([0 Inf 10^-7 Inf])
    xlabel('P');
    ylabel('reconstructieerror');
    print([file_for 'Elog' '.pdf'],'-dpdf')

    figure;
    errorbar(P,mp,sp);
    hold on;
    plot(P,P,'k--');
    hold off;
    xlabel('gegeven P');
    ylabel('ontmegde P');
    axis([0 Inf 0 Inf])
    axes('position',[.2 .60 .25 .25])
    hold on;
    errorbar(smallP,mp(1:zoom_points),sp(1:zoom_points));
    plot(smallP,smallP,'k--');
    hold off;
    axis([0 Inf 0 Inf])
    print([file_for 'PP' '.pdf'],'-dpdf')
    
    figure;
    errorbar(P,ma,sa);
    xlabel('P');
    ylabel('hoek');
    axis([0 Inf 0 Inf])
    axes('position',[.60 .60 .25 .25])
    errorbar(smallP,ma(1:zoom_points),sa(1:zoom_points));
    axis([0 Inf 0 Inf])
    print([file_for 'angle' '.pdf'],'-dpdf')
    
    figure;
    h1 = errorbar(mp,me,se);
    hold on;
    h2 = errorbar(mp,ml,sl,'r');
    hold off;
    set(get(h1,'Parent'), 'YScale', 'log')
    set(get(h2,'Parent'), 'YScale', 'log')
    axis([0 Inf 10^-7 Inf])
    xlabel('gemeten P');
    ylabel('reconstructieerror');
    print([file_for 'EPlog' '.pdf'],'-dpdf')
   
end

function [mcorr,scorr,merr,serr,mrefl,srefl,mangle,sangle,mlin,slin] = mean_correctness(P,hits)
    [corr,err,refl,angle,lin] = linfixfast(repmat(P(:),1,hits));
    [mcorr,scorr] = grow_mats(corr,P);
    [merr,serr] = grow_mats(err,P);
    [mrefl,srefl] = grow_mats(refl,P);
    [mangle,sangle] = grow_mats(angle,P);
    [mlin,slin] = grow_mats(lin,P);
end

function [m,s] = grow_mats(mat,P)
m = reshape(mean(mat.'),size(P));
s = reshape(std(mat.'),size(P));
end

function [tab,corr,error] = linfix(P)
corr = P*0;
error = P*0;
[library,~,~] = load_Alina();
for i=1:size(P(:),1)
    [pixels,solution] = mult_test_alina(library,P(i));

    [indices, abundance, rP, ~, ~, ~, err] = libraryperm(pixels,library,true,false,true,false);
    [ uma, umP, ~, ~, ume ] = multilin_psmall( pixels, solution.endmem,true,false);

    disp(indices);
    disp(solution.indices);
    tab{i} = table([abundance;uma;solution.a],[rP;umP;solution.P],[indices;solution.indices;solution.indices],[err;ume;0],...
        'RowNames',{'Selected','Unmixed','Exact'});
    corr(i) = mean((indices(:)-solution.indices(:))==0);
    error(i) = err;
end



end

function [corr,error,refl,angle,linerror] = linfixfast(P)
pixels = zeros(size(P(:),1),53);
compare = zeros(size(P(:),1),4);
[library,~,~,rl] = load_Alina();
for i=1:size(P(:),1)
    [pixels(i,:),solution] = mult_test_alina(library,P(i));
    compare(i,:) = solution.indices;
end
%[indices, abundance, rP, ~, ~, ~, error] = libraryperm(pixels,library,true,false,true,false);
[indices,error,linerror,refl] = library_fast(library,rl,pixels);
%[ uma, umP, ~, ~, ume ] = multilin_psmall( pixels, solution.endmem,true,false);

disp(indices);
disp(compare.');
disp((indices-compare.')==0);
%tab{i} = table([abundance;uma;solution.a],[rP;umP;solution.P],[indices;solution.indices;solution.indices],[err;ume;0],...
%    'RowNames',{'Selected','Unmixed','Exact'});
corr = mean(((indices-compare.')==0));

corr = reshape(corr,size(P));
error = reshape(error,size(P));
linerror = reshape(linerror,size(P));
refl = reshape(refl,size(P));
angle = reshape(max(get_angle(indices,compare.',library)),size(P));

end

function angle = get_angle(indices,compare,library)
angle = indices*0;
[libcount,endcount] = size(indices);
for lib=1:libcount
    book = library{lib};
    for endm=1:endcount
        real = compare(lib,endm);
        mesured = indices(lib,endm);
        if real==0 || mesured == 0 || real==mesured
            angle(lib,endm) = 0;
        else
            realspec = book(real,:);
            mesuredspec = book(mesured,:);
            angle(lib,endm) = acos(dot(realspec,mesuredspec)/(norm(realspec)*norm(mesuredspec)));
        end
    end
end


end

function [table,errors] = table_show_result()
table = {};
for i=1:20
    fieldname = ['tablerow' num2str(i)];
    loadtable = load('tabletestS',fieldname);
    if isfield(loadtable,fieldname)
        table{i} = loadtable.(fieldname);
    end
end

m = length(table);
errors = zeros(5,m);
runtime = zeros(5,m);

for i=1:m
    errors(:,i) = table{i}.error(:);
    runtime(:,i) = table{i}.runtime(:);
end
close all;
bar(1:m,errors.')
set(gca, 'YScale', 'log')
legend({'Linair';'Semiliniair (P)';'Semiliniair (Pi)';'Multilineair (P)';'Multilineair (Pi)'})
xlabel positie
ylabel reconstructieerror
print('exp_errors','-dpdf')
figure;
bar(1:m,runtime.')
set(gca, 'YScale', 'log')
legend({'Linair';'Semiliniair (P)';'Semiliniair (Pi)';'Multilineair (P)';'Multilineair (Pi)'})
xlabel positie
ylabel runtime
print('exp_runtime','-dpdf')
figure;
bar(1:m,runtime.')
legend({'Linair';'Semiliniair (P)';'Semiliniair (Pi)';'Multilineair (P)';'Multilineair (Pi)'})
xlabel positie
ylabel runtime
print('exp_runtime_honest','-dpdf')
end

function table = Table_test_Alina(row)

save('tabletestS')

for i=1:21
    table.(['tablerow' num2str(i)]) = pixelcompare(i,row);
    save('tabletestS','-Struct','table','-append')
    
end

end


function p = pixelcompare(x,y)
[libs,~,data] = load_Alina();
pixel = data(y,x,:)
abundance = zeros(5,4);
reflections = zeros(5,4);
indices = zeros(5,4);
P = zeros(5,1);
error = zeros(5,1);
runtime = zeros(5,1);
tic()
[indices(1,:),abundance(1,:), P(1), reflections(1,:), ~, ~, error(1)] = libraryperm(pixel,libs,true,true,false,false);
runtime(1) = toc();
tic()
[indices(2,:),abundance(2,:), P(2), reflections(2,:), ~, ~, error(2)] = libraryperm(pixel,libs,true,false,false,false);
runtime(2) = toc();
tic()
[indices(3,:),abundance(3,:), P(3), reflections(3,:), ~, ~, error(3)] = libraryperm(pixel,libs,true,false,true,false);
runtime(3) = toc();
tic()
[indices(4,:),abundance(4,:), P(4), reflections(4,:), ~, ~, error(4)] = libraryperm(pixel,libs,false,false,false,false);
runtime(4) = toc();
tic()
[indices(5,:),abundance(5,:), P(5), reflections(5,:), ~, ~, error(5)] = libraryperm(pixel,libs,false,false,true,false);
runtime(5) = toc();
RowName = {'Linair';'Mixed (P)';'Mixed (Pi)';'Multilineair (P)';'Multilineair (Pi)'};
p = table(runtime,error,abundance,P,reflections,indices,'RowNames',RowName);

end

%function reccompare();



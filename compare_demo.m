function c = compare_demo()
c.testpix = @testpix;
c.show = @show_table;
c.c = @cdata;
c.cf = @cfunc;
c.lfpd = @save_latex_file_pixel_data;
end


function [tabl,error,runtime,abundance,reflectance,names] = testpix(x,y)

[libs,~,data,rl] = load_Alina();
pixel = reshape(data(y,x,:),[1,53]);


% allsystems = { ...
%     multilin_options(1,true,false,false,600), ...
%     multilin_options(1,false,false,false,600), ...
%     multilin_options(1,false,true,false,600), ...
%     multilin_options(1,false,false,true,600), ...
%     multilin_options(1,false,true,true,600), ...
%     };

allsystems = { ...
    %multilin_options(0,false,false,false,600), ...
    %multilin_options(0,false,false,true,600), ...
    %multilin_options(0,false,true,false,600), ...
    %multilin_options(0,false,true,true,600), ...
    multilin_options(1,true,false,false,600), ...
    multilin_options(1,false,true,false,600), ...
    multilin_options(1,false,true,true,600), ...
    multilin_options(1,false,false,true,600), ...    
    multilin_options(1,false,false,false,600), ...
    multilin_options(2,true,false,false,600), ...
    multilin_options(2,false,true,false,600), ...
    multilin_options(2,false,true,true,600), ...
    multilin_options(2,false,false,true,600), ...    
    multilin_options(2,false,false,false,600), ...
    };


    
l = length(allsystems);
abundance = zeros(l,4);
reflectance = zeros(l,4);
error = zeros(l,1);
runtime = zeros(l,1);
freedom = zeros(l,1);
names = {};
markers = cell(l,1);

i = 0;
for op  = allsystems
    i = i + 1;
    options = op{1};
    disp(options)
    names{i} = options.name;
    markers{i,1} = options.marker;
    
    tic()

    switch options.select
        case 0
            indices = library_custom(pixel,libs,options);
        case 1
            indices = MESMA_brute_small(pixel.',rl);
        case 2
            addpath('AAM/')
            indices = AAM(pixel.',rl);
        otherwise
            assert(true,'Invalid method')
    end
    
    perm = find(indices);
    endmembers = build_endmem(libs,indices);

    if options.linsolve
        size(abundance(i,:))
        [a,~,error(i)] = FCLSU(endmembers,pixel);
        abundance(i,perm) = a.';
    else
        [ abundance(i,perm), ~, reflectance(i,perm), ~, error(i) ] = multilin_custom(pixel, endmembers,options);

    end
    runtime(i) = toc();
    freedom(i,:) = options.dof;
end

tabl = table(runtime,error,freedom,abundance,reflectance,markers,'RowNames',names);

end

function [tabl,runtime,error] = testarea(x,y,w,h)
error = [];
runtime = [];

for i=1:w
    for j=1:h
        [t,e,r] = testpix(x+i-1,y+j-1);
        error = [error , e];
        runtime = [runtime , r];
    end
end
names = t.Properties.RowNames;

error = mean(error.').'
runtime = mean(runtime.').'
tabl = table(runtime,error,'RowNames',names)

end


function show_table(numberofsystems,rownumber)
errors = zeros(numberofsystems,19);
runtime = errors;


for x=1:19
    [~,errors(:,x),runtime(:,x),~,~,names] = testpix(x,rownumber);
end


close all;
bar(1:19,errors.')
set(gca, 'YScale', 'log')
legend(names)
xlabel positie
ylabel reconstructieerror
colormap jet
print('full_errors','-dpdf')
figure;
bar(1:19,runtime.')
set(gca, 'YScale', 'log')
legend(names)
xlabel positie
ylabel reconstructieerror
colormap jet
print('full_runtimes','-dpdf')

end

function endmem = build_endmem(lib,indices)
perm = find(indices);
len = size(perm(:),1);
bands = size(lib{1},2);
endmem = zeros(len,bands);
for i=1:len
    endmem(i,:) = lib{perm(i)}(indices(perm(i)),:);
end

end
function cdata()

trees = testarea(1,3,3,1);
mix1 = testarea(4,3,3,1);
road = testarea(8,3,3,1);
mix2 = testarea(11,3,3,1);
curb = testarea(13,3,3,1);
grass = testarea(16,3,3,1);

names = trees.Properties.RowNames;

errors = table(trees.error,mix1.error,road.error,mix2.error,curb.error,grass.error,'RowNames',names);
errorsrn = table(rne(trees),rne(mix1),rne(road),rne(mix2),rne(curb),rne(grass),'RowNames',names);
runtime = table(trees.runtime,mix1.runtime,road.runtime,mix2.runtime,curb.runtime,grass.runtime,'RowNames',names);

save_latex_file(errors,'TabErrorsA');
save_latex_file(errorsrn,'TabErrorsrnA');
save_latex_file(runtime,'TabRuntimeA');
%save_latex_file(table2,'TabRoad');
%save_latex_file(table3,'TabGrass');
save('meantables.mat','errors','errorsrn','runtime');



%show_table(11,3);
%cfunc()
end

function rn = rne(tab)

rn = tab.error/max(tab.error);

end

function cfunc()
load('fulltables.mat')
figure;
table_grapher(trees)
title('Trees')
print('plot_trees','-dpdf')
figure;
table_grapher(road)
title('Road')
print('plot_road','-dpdf')
figure;
table_grapher(grass)
title('Grass')
print('plot_grass','-dpdf')
figure;
for i=1:size(trees.runtime)
    plot(i,trees.freedom(i),trees.markers{i})
    hold on;
end
ylabel('degrees of freedom')

set(gca,'XTick',[])
set(gca,'XColor','w')
lh = legend(trees.Properties.RowNames);
set(lh,'location','northeastoutside');
print('plot_freedom','-dpdf')
figure;
for i=1:size(trees.runtime)
    plot(1,1,trees.markers{i})
    hold on;
end

axis off;
lh = legend(trees.Properties.RowNames);
set(lh,'location','northeastoutside');
print('plot_legend','-dpdf')
end





%Deze code snipped komt uit de example file van latexTable
function save_latex_file_pixel_data(table,filename)
table.markers = [];
input.data = table;
input.dataFormat = {'%.3f',2,'%.0f',1,'%.3f',8};
input.tableColLabels = {'runtime','error','freedom','$a_1$','$a_2$','$a_3$','$a_4$','$P_1$','$P_2$','$P_3$','$P_4$'};
latex = latexTable(input);
% save LaTex code as file


fid=fopen( [filename '.tex' ],'w');
[nrows,~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
fprintf('\n... your LaTex code has been saved in your working directory\n');
   
end
%Deze code snipped komt uit de example file van latexTable
function save_latex_file(table,filename)
input.data = table;
latex = latexTable(input);
% save LaTex code as file
fid=fopen( [filename '.tex' ],'w');
[nrows,~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
fprintf('\n... your LaTex code has been saved in your working directory\n');
   
end


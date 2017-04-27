function c = compare_demo()
c.testpix = @testpix;
c.show = @show_table;
c.c = @cdata;
c.cf = @cfunc;
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
    multilin_options(0,false,false,false,600), ...
    multilin_options(0,false,false,true,600), ...
    multilin_options(0,false,true,false,600), ...
    multilin_options(0,false,true,true,600), ...
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

tabl = table(runtime,error,freedom,abundance,reflectance,markers,'RowNames',names)

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

trees = testpix(3,3);
road = testpix(11,3);
grass = testpix(17,3);

%save_latex_file(table1,'TabTree');
%save_latex_file(table2,'TabRoad');
%save_latex_file(table3,'TabGrass');
save('fulltables.mat','trees','road','grass');

%show_table(11,3);
cfunc()
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


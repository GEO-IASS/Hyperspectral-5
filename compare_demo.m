function c = compare_demo()
c.testpix = @testpix;


end


function [tabl,error,runtime,abundance,reflectance] = testpix(x,y)

[libs,~,data,rl] = load_Alina();
pixel = reshape(data(y,x,:),[1,53]);


allsystems = { ...
    %multilin_options(false,false,false,false,600), ...
    multilin_options(true,true,false,false,600), ...
    multilin_options(true,false,false,false,600), ...
    multilin_options(true,false,true,false,600), ...
    multilin_options(true,false,false,true,600), ...
    multilin_options(true,false,true,true,600)};
    
l = length(allsystems);
abundance = zeros(l,4);
reflectance = zeros(l,4);
error = zeros(l,1);
runtime = zeros(l,1);
freedom = zeros(l,1);
names = {};


i = 0;
for op  = allsystems
    i = i + 1;
    options = op{1};
    disp(options)
    names{i} = options.name;
    tic()

    if options.linapprox
        indices = library_fast(libs,rl,pixel);
    else
        indices = library_custom(pixel,libs,options);
    end
    
    perm = find(indices);
    endmembers = build_endmem(libs,indices);

    if options.linsolve
        size(abundance(i,:))
        [a,~,error(i)] = SCLSU(endmembers,pixel);
        abundance(i,perm) = a.';
    else
        [ abundance(i,perm), ~, reflectance(i,perm), ~, error(i) ] = multilin_custom(pixel, endmembers,options);

    end
    runtime(i) = toc();
    freedom(i,:) = options.dof;
end

tabl = table(runtime,error,freedom,abundance,reflectance,'RowNames',names)

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
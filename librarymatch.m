function [sub_indices, abundance, P, reflectance, ...
    reconstruct, albedos, error] = librarymatch(pixels,libs,cutoff,linapprox,independant_P,shadow)
    len = length(libs);
    lsize = zeros(1,len);
    for i=1:len
        lsize(i) = size(libs{i},1);
        bands = size(libs{i},2);
    end
    min(size(libs,2));
    counter = lsize;
    sub_indices = counter*0;
    error = cutoff;
    cumsize = cumprod(lsize)./lsize;
    numberofcalculations = sum((counter-1).*cumsize);
    %disp(['confs left: ' num2str(numberofcalculations)])
    tic_disp = tic();
    tic_loop = tic();
    while true
        if linapprox
            [~,~,err] = SCLSU(build_endmem(libs,counter),pixels);
        else
            [~,~,~,~,err] = multilin_psmall(pixels,build_endmem(libs,counter),independant_P,shadow);
        end
        
        
        %disp(['current error: ' num2str(err)]);
        if err < error
            sub_indices = counter;
            error = err;
        end
        %disp([num2str(counter) '->'  num2str(err)])
        nonzeros = find(counter-1);
        if numel(nonzeros) == 0
            break
        end
        f = nonzeros(1);
        if toc(tic_disp) > 5
            calcleft = sum((counter-1).*cumsize);
            timeleft = toc(tic_loop)*calcleft/(numberofcalculations-calcleft);
            disp(['Time elapsed: ' datestr(datenum(0,0,0,0,0,toc(tic_loop)),'dd:HH:MM:SS')])
            disp(['Time left:~ ' datestr(datenum(0,0,0,0,0,timeleft),'dd:HH:MM:SS')])
            disp(['confs left: ' num2str(calcleft)])
            disp(['current error: ' num2str(error)])
            tic_disp = tic();
        end
        if f > 1
            
            for i=1:f-1
                counter(i) = lsize(i);
            end
        end
        counter(f) = counter(f) - 1;
        
    end 
    if sum(sub_indices) == 0
        %disp('no valid configurations')
        sub_indices = sub_indices + 1;
    else
        %disp(['found configuration: ' num2str(sub_indices)])
        %disp(['configuration error: ' num2str(error)])
    end
    %disp('performing multilineair unmixing')
    %[abundance, P, reflectance, reconstruct, albedos, error] = ...
    %        multilin(pixels,build_endmem(libs,sub_indices),Inf(),Inf());
    [abundance,reconstruct,error] = SCLSU(build_endmem(libs,sub_indices),pixels);
    abundance = abundance.';
    P = 0;
    reflectance = 0;
    albedos = 0;
    %disp(['ERROR: ' num2str(error)])
end


function endmem = build_endmem(lib,indices)
len = length(lib);
bands = size(lib{1},2);
endmem = zeros(len,bands);
for i=1:len
    endmem(i,:) = lib{i}(indices(i),:);
end

end
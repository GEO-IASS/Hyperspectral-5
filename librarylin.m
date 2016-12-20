function [lib_index, sub_indices, abundance, P, reflectance, ...
    reconstruct, albedos, error] = librarylin(pixels,libs,maxmat)
tic()
len = length(libs);
errors = ones(1,len);

workers = 4
if isempty(gcp('nocreate'))
    parpool(workers)
end


disp('error   size*error sub_indices')
for l=1:len
    disp(['    LIBRARY ' num2str(l)])
    biglib = libs{l};
    [count,~] =size(biglib);

    [permutes{l}, errors(l)] = unpack_lib_hybrid(biglib,pixels,maxmat,64);

end




[~, lib_index] = min(errors);
sub_indices = permutes{lib_index};
delete(gcp)
lib = libs{lib_index};
smlib = lib(sub_indices,:);
[abundance, P, reflectance, reconstruct, albedos, error] = ...
            multilin(pixels,smlib);
disp(['Time elapsed: ' datestr(datenum(0,0,0,0,0,toc()),'HH:MM:SS')])

end

function [perm,error] = unpack_lib_serial(biglib,pixels,maxsize)
    
    errors = Inf(1,2^size(biglib,1)-1);
    [count,~] =size(biglib);
    [ab, ~, ~, ~, ~, err] = multilin(pixels,biglib);
    for p=1:2^count
        errors(p) = check_perm(biglib,p,count,maxsize,pixels);
    end

    [error, p] = min(errors);
    perm = find(dec2bin(p,count) - '0');
end

function [perm,error] = unpack_lib_parallel(biglib,pixels,maxsize)
    
    errors = Inf(1,2^size(biglib,1)-1);
    [count,~] =size(biglib);
    [ab, ~, ~, ~, ~, err] = multilin(pixels,biglib);
    parfor p=1:2^count
        errors(p) = check_perm(biglib,p,count,maxsize,pixels);
    end

    [error, p] = min(errors);
    perm = find(dec2bin(p,count) - '0');
end


function [perm,error] = unpack_lib_hybrid(biglib,pixels,maxsize,batch)
    errors = Inf(1,2^size(biglib,1)-1);
    [count,~] =size(biglib);
    for ps=1:batch:2^count
        disp(['starting batch: ' num2str(ps) '/' num2str(2^count)])
        parfor p=ps:ps+batch-1
            errors(p) = check_perm(biglib,p,count,maxsize,pixels);
        end
    end

    [error, p] = min(errors);
    perm = find(dec2bin(p,count) - '0');
end



function [error, perm] = check_perm(biglib,index,count,maxsize,pixels)
    perm = find(dec2bin(index,count) - '0');
    psize = size(perm,2);
    if psize <= 1 || psize > maxsize || index>=2^count
        error = Inf;
        return
    end
    smlib = biglib(perm,:);
    [~, ~, ~, ~, ~, err] = multilin(pixels,smlib,0.05);
    disp([num2str(2^count-index, '%04.0f') '-> ' num2str(err, '%07.4f') ', '  ...
        num2str(err*psize, '%07.4f') ': ' num2str(perm)  ])
    error = err*psize;

end
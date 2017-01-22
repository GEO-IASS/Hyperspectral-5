function [indices,error] = library_custom(pixel,libs,options)

len = length(libs);
error = Inf();

for p=1:2^len-1
    perm = find(dec2bin(p,len) - '0');
    disp(['using libraries ' num2str(perm)])
    slib = libs(perm);
    len = length(slib);
    lsize = zeros(1,len);
    for i=1:len
        lsize(i) = size(slib{i},1);
    end
    counter = lsize;
    while true
        [~,~,~,~,err] = multilin_custom(pixel,build_endmem(slib,counter),options);
        nonzeros = find(counter-1);
                
        if err < error
            sub_indices = counter;
            indices = indices*0;
            indices(perm) = sub_indices;
            error = err;
        end
        if numel(nonzeros) == 0
            break
        end
        f = nonzeros(1);
        if f > 1
            
            for i=1:f-1
                counter(i) = lsize(i);
            end
        end
        counter(f) = counter(f) - 1;
    end
end

end


function endmem = build_endmem(lib,indices)
len = length(lib);
bands = size(lib{1},2);
endmem = zeros(len,bands);
for i=1:len
    endmem(i,:) = lib{i}(indices(i),:);
end

end
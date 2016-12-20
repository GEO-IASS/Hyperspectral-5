function [BIGindices, abundance, P, reflectance, ...
    reconstruct, albedos, error] = libraryperm(pixels,libs,linapprox,linmix,independant_P,shadow)
    len = length(libs);
    lsize = zeros(1,len);
    maxmat = Inf();
    maxerror = Inf();
    
    [pixcount,bands] = size(pixels);
    abundance = zeros(pixcount,len);
    P = zeros(pixcount,1);
    reconstruct = zeros(pixcount,bands);
    error = zeros(1,pixcount);
    BIGindices = zeros(pixcount,len);
    for pixelindex=1:pixcount
        sub_indices = zeros(1,len);
        disp(['pixels: ' num2str(pixelindex) '/' num2str(pixcount)])
        pixel = pixels(pixelindex,:);
        cutoff = maxerror;
        for p=1:2^len-1
            perm = find(dec2bin(p,len) - '0');
            disp(['using libraries ' num2str(perm)])
            psize = size(perm,2);
            if psize > maxmat
                continue
            end
            slib = libs(perm);
            [indices, ~, ~, ~, ~, ~,err] = librarymatch(pixel,slib,Inf(),linapprox,independant_P,shadow);
            %disp(['error: ' num2str(err)])
            if cutoff > err
                sub_indices = zeros(1,len);
                sub_indices(perm) = indices;
                cutoff = err;
            end
        end
        %disp(['ERROR: ' num2str(cutoff)])
        %disp(['indices: ' num2str(sub_indices)])
        perm = find(sub_indices);
        if linmix
            [pix_abundance,pix_reconstruct,pix_error] = SCLSU(build_endmem(libs,sub_indices,perm),pixel);
            pix_P = 0;
        else
            [pix_abundance,pix_P,~,pix_reconstruct,pix_error] = multilin_psmall(pixel,build_endmem(libs,sub_indices,perm),independant_P,shadow);
        end
        
        local_abundance = zeros(1,len);
        local_abundance(perm) = pix_abundance.';
        abundance(pixelindex,:) = local_abundance;
        %abundance(pixelindex,perm) = pix_abundance.';
        P(pixelindex,:) = pix_P;
        reconstruct(pixelindex,:) = pix_reconstruct.';
        BIGindices(pixelindex,:) = sub_indices;
        error(pixelindex) = pix_error;
    end
    albedos=0;
    reflectance=0;

    
   
end

function endmem = build_endmem(lib,indices,perm)
len = size(perm(:),1);
bands = size(lib{1},2);
endmem = zeros(len,bands);
for i=1:len
    endmem(i,:) = lib{perm(i)}(indices(perm(i)),:);
end

end
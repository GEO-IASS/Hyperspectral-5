function [abundance,reconstruct,err] = SCLSU(endmembers,spectrum)
speccol = spectrum.';
[bands,pixels] = size(speccol);
[endc,bands] = size(endmembers);
abundance = zeros(endc,pixels);
if endc~=1
    shadow = endmembers(1,:);
    endmembers_noshad = endmembers(2:end,:) - repmat(shadow,endc-1,1);


    abundance(2:end,:) = pinv(endmembers_noshad.')* (speccol - repmat(shadow,pixels,1).');
    
    abundance(1,:) = 1-sum(abundance(2:end,:));
else
    abundance = ones(1,pixels);
end
reconstruct = endmembers.'*abundance;

%check for nonnegativity
if all(abundance > 0)
    err = sumsqr(speccol(:)-reconstruct(:));
else
    err = Inf();
end
end

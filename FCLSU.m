function [abundance,reconstruct,err] = FCLSU(endmembers,spectrum)
speccol = spectrum.';
abundance = zeros(size(endmembers.',2),size(speccol,2));
[bands,pixels] = size(speccol);
[endc,bands] = size(endmembers);
for i=1:size(speccol,2)
    shadow = endmembers(1,:);
    endmembers_noshad = endmembers(2:end,:) - repmat(shadow,endc-1,1);


    abundance(2:end,:) = lsqnonneg(endmembers_noshad.',speccol - repmat(shadow,pixels,1).');
    
    abundance(1,:) = 1-sum(abundance(2:end,:));    
    
reconstruct = endmembers.'*abundance;

err = sumsqr(speccol(:)-reconstruct(:));
end

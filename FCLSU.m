function [abundance,reconstruct,err] = FCLSU(endmembers,spectrum)
speccol = spectrum.';
abundance = zeros(size(endmembers.',2),size(speccol,2));
for i=1:size(speccol,2)
    abundance(:,i) = lsqnonneg(endmembers.',speccol(:,i));
    
    
reconstruct = endmembers.'*abundance;

err = sumsqr(speccol(:)-reconstruct(:));
end

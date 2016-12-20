function [pixels, endmembers,solution] = build_test_spectrum(pixelcount,endmembers_exact,noise,P)
[emcount,bandcount] = size(endmembers_exact);
%albedos = rand(emcount,bandcount);

%memrefl = rand(emcount,1);
%memrefl_exp = repmat(memrefl,1,bandcount);
%endmembers_exact = ((1 - memrefl_exp).* albedos)./(1 - memrefl_exp.* albedos);
endmembers = normrnd(endmembers_exact,noise);


pixrefl = P;
abundances = randa(emcount*pixelcount,'double');
abundances = reshape(abundances,[emcount pixelcount]);

%abundances = bsxfun(@rdivide, abundances, sum(abundances)); %normalises the abundances
%solution.R = memrefl.';
solution.P = pixrefl.';
solution.a = abundances.';
solution.endmem = endmembers_exact;
%solution.w = albedos;
pixels_exact = zeros(pixelcount,bandcount);
for i=1:pixelcount
    abundances_exp = repmat(abundances(:,i),1,bandcount);
    pixels_exact(i,:) = (1 - pixrefl)*sum(abundances_exp.* endmembers_exact)./(1 - pixrefl*sum(abundances_exp.* endmembers_exact));
solution.pix = pixels_exact;
pixels = normrnd(pixels_exact,noise);

end
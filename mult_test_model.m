function [pixels, library,solution] = mult_test_model(emcount,libsize,bandcount,P)


library = generate_libs(create_spectrum(emcount,bandcount,5),libsize);
lib_indices = randi([1 libsize],1,emcount);
%plot(1:bandcount,libs{1})
endmembers = zeros(emcount,bandcount);
for i=1:emcount
    loclib = library{i};
    endmembers(i,:) = loclib(lib_indices(i),:);
end
[pixels, ~, solution] = build_test_spectrum(1,endmembers,0,P);
solution.indices = lib_indices;

end

function spec = create_spectrum(speccount,bandcount,variance)

overhead = rand(speccount,bandcount + variance - 1);
spec = zeros(speccount,bandcount);


for i=1:variance
    spec = spec + overhead(:,i:i+bandcount-1);
end
spec = spec/variance;
    
end

function libs = generate_libs(spectra,libsize,spike)

libs = {};

for i=1:size(spectra,1)
    libs{i} = binornd(spike,repmat(spectra(i,:),libsize,1))/spike;
end

end
function [pixels,solution] = mult_test_alina(library,P)


%library = generate_libs(create_spectrum(emcount,bandcount,5),libsize,100);

emcount = length(library);
bandcount = size(library{1},2);
libsize = zeros(1,emcount);
for i=1:emcount
    libsize(i) = size(library{i},1);
end
    
lib_indices = floor(rand(size(libsize)) .* libsize)+1;
%plot(1:bandcount,libs{1})
endmembers = zeros(emcount,bandcount);
for i=1:emcount
    loclib = library{i};
    endmembers(i,:) = loclib(lib_indices(i),:);
end
[pixels, ~, solution] = build_test_spectrum(1,endmembers,0,P);
solution.indices = lib_indices;

end


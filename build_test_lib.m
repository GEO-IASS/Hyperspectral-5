function [pixels, library,solution] = build_test_lib(pixelcount,bandcount,libcount,libsize,indices,noise)
[pixels, endmem ,solution] = build_test_spectrum(pixelcount,bandcount,size(indices,2),noise);
for i=1:libcount
    library{i} = rand(libsize,bandcount);

    library{i}(indices(i),:) = endmem(i,:);


end

end
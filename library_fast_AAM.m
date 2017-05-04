function [indices,error,linerror,P] = library_fast_AAM(libs,rl,pixels)
disp('Selecting')
[indices,~,~,linerror] = AAM(pixels.',rl,1);
linerror = linerror.^2;
pixcount = size(indices,2);
error = zeros(1,pixcount);
P = zeros(1,pixcount);
disp('unmixing')
for i=1:pixcount
    if mod(i,1000) == 0
        disp(i/pixcount) 
    end
    [~,P(i),~,~,error(i)] = multilin_psmall(pixels(i,:),build_endmem(libs,indices(:,i)),false,false);
end


end

function endmem = build_endmem(lib,indices)
perm = find(indices);
len = size(perm(:),1);
bands = size(lib{1},2);
endmem = zeros(len,bands);
for i=1:len
    endmem(i,:) = lib{perm(i)}(indices(perm(i)),:);
end

end
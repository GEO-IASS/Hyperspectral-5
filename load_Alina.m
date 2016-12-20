function [lib,data,Xim,rawlib] = load_Alina()
load('Alina/demo.mat');
rawlib = endmembers;
lib = cellfun(@transpose,endmembers,'UniformOutput',false);
[x,y,b] = size(Xim);
data = reshape(Xim,x*y,b);

end
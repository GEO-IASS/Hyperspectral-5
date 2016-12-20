function toterr = multilin_errorguess( pixels, endmembers,cutoff,tolerance)
%MULTILIN Summary of this function goes here
%   Detailed explanation goes here
[endcount,~] = size(endmembers);
[pixelcount,bands] = size(pixels);

default = zeros(2,endcount);
default(1,:) = ones(1,endcount)/endcount;
default(2,:) = ones(1,endcount)/endcount/3;

toterr = 0;
for i=1:pixelcount
    spectrum = pixels(i,:);
    [inputmat,fval] = fmincon(@(I) guassian_error_lim(I,spectrum,endmembers), ...
                default,[],[],[],[],zeros(2,endcount),ones(2,endcount),[], ...
                optimset('MaxFunEvals',600, 'Display','off', ...
                'algorithm','interior-point','TolFun',tolerance));
    toterr = toterr + fval/pixelcount;
    if toterr > cutoff
        return 
    end
end
end

function albedos = get_albedos(R,endmembers)
[endcount,bands] = size(endmembers);
exp_R = repmat(R.',1,bands);

albedos = endmembers./(endmembers.*exp_R + 1 - exp_R);


end

function spectrum = get_spectrum(b,c,albedos)

b_sum = b * albedos;
c_sum = c * albedos;

spectrum = b_sum./(1-c_sum);
end

function [error] = guassian_error_lim(inputmat, spectrum, endmembers)
b = inputmat(1,:);
c = inputmat(2,:);

% inefficient version
% albedos = zeros([endcount,bands]);
% for i=1:endcount
%    albedos(i,:) = endmembers(i,:)./( R(i) * endmembers(i,:) + (1 - R(i)));
% end

%efficient version
albedos = get_albedos(zeros(size(b)),endmembers);

build_spectrum = get_spectrum(b,c,albedos);

error = norm(spectrum - build_spectrum,2);
end
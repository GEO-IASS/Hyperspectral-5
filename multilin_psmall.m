function [ abundance, P, reflectance, reconstruct, error ] = multilin_psmall( pixels, endmembers,independant_P,shadow)
%MULTILIN Summary of this function goes here
%   Detailed explanation goes here
%   dependant_P says if the P value should be endmember dependant
%   allowshadow says if the 
[endcount,~] = size(endmembers);
[pixelcount,bands] = size(pixels);

%sunlight = mean(pixels)
%endmembers = endmembers./repmat(sunlight,size(endmembers,1),1);
%<- werkt niet en is onnodig

abundance = zeros(pixelcount,endcount);
reflectance = zeros(pixelcount,endcount);
P = zeros(pixelcount,1);
reconstruct = pixels;
errors = zeros(pixelcount,1);
shadows = zeros(pixelcount,1);
default = zeros(2,endcount);


toterr = 0;
R = zeros(1,endcount);


for i=1:pixelcount
    spectrum = pixels(i,:);
    default(1,:) = SCLSU(endmembers,pixels(i,:));
    default(2,:) = ones(1,endcount)*0;
    inputmat = fmincon(@(I) guassian_error_lim(I,R,spectrum,endmembers,independant_P,shadow), ...
                default(:,:),[],[],[],[],zeros(2,endcount),ones(2,endcount),[], ...
                optimset('MaxFunEvals',600, 'Display','off', ...
                'algorithm','interior-point'));
    b = inputmat(1,:);
    c = inputmat(2,:);
    %R = zeros(1,endcount)
    
    if ~shadow
        %sum to one constraint
        totalsum = sum(b) + sum(c);
        b = b/totalsum;
        c = c/totalsum;
    end
    
    %independant P constrained
    if ~independant_P
        a_ = (b+c);
        P_ = sum(c);
        b = a_*(1-P_);
        c = a_*P_;
    end
    

    shadows(i) = (1-sum(c))/sum(b);
    abundance(i,:) = shadows(i)*b+c;
    reflectance(i,:) = c./abundance(i,:);
    P(i) = sum(c);
    reconstruct(i,:) = get_spectrum(b,c,get_albedos(R,endmembers));

    errors(i) = sumsqr(spectrum-reconstruct(i,:));
end
error = sum(errors);

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

function [error] = guassian_error_free(inputmat, spectrum, endmembers)

b = inputmat(1,:);
c = inputmat(2,:);
R = inputmat(3,:);

% inefficient version
% albedos = zeros([endcount,bands]);
% for i=1:endcount
%    albedos(i,:) = endmembers(i,:)./( R(i) * endmembers(i,:) + (1 - R(i)));
% end

%efficient version
albedos = get_albedos(R,endmembers);

reconstruct = get_spectrum(b,c,albedos);

error = sumsqr(spectrum - reconstruct);

end


function [error] = guassian_error_lim(inputmat, R, spectrum, endmembers,independant_P,shadow)
b = inputmat(1,:);
c = inputmat(2,:);

if ~shadow
    %sum to one constraint
    totalsum = sum(b) + sum(c);
    b = b/totalsum;
    c = c/totalsum;
end

%independant P constrained
if ~independant_P
    a_ = (b+c);
    P_ = sum(c);
    b = a_*(1-P_);
    c = a_*P_;
end

% inefficient version
% albedos = zeros([endcount,bands]);
% for i=1:endcount
%    albedos(i,:) = endmembers(i,:)./( R(i) * endmembers(i,:) + (1 - R(i)));
% end

%efficient version
albedos = get_albedos(R,endmembers);

reconstruct = get_spectrum(b,c,albedos);

error = sumsqr(spectrum - reconstruct);
end
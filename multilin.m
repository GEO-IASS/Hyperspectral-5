function [ abundance, P,reflectance, reconstruct, albedos, error ] = multilin( pixels, endmembers,cutoff, linfac)
%MULTILIN Summary of this function goes here
%   Detailed explanation goes here
[endcount,~] = size(endmembers);
[pixelcount,bands] = size(pixels);

%sunlight = mean(pixels)
%endmembers = endmembers./repmat(sunlight,size(endmembers,1),1);
%<- werkt niet en is onnodig

abundance = zeros(pixelcount,endcount);
reflectance = zeros(pixelcount,endcount);
P = zeros(pixelcount,1);
reconstruct = pixels;
albedos = endmembers*0;
errors = zeros(pixelcount,1);
shadows = zeros(pixelcount,1);
Rstruct = zeros(pixelcount,endcount);

default = zeros(3,endcount);
default(1,:) = ones(1,endcount)/endcount;
default(2,:) = ones(1,endcount)/endcount/3;
default(3,:) = ones(1,endcount)/endcount/8;
guesserr = SCLSU(endmembers,pixels)/pixelcount;
if guesserr > cutoff*linfac
    error = Inf();
    return 
end
toterr = 0;
for i=1:pixelcount
    spectrum = pixels(i,:);
    [inputmat,fval] = fmincon(@(I) guassian_error_free(I,spectrum,endmembers), ...
                default,[],[],[],[],zeros(3,endcount),ones(3,endcount),[], ...
                optimset('MaxFunEvals',6000, 'Display','off', ...
                'algorithm','interior-point'));
    toterr = toterr + fval;
    if toterr > cutoff
        error = Inf();
        return 
    end
    %b = inputmat(1,:);
    %c = inputmat(2,:);
    R = inputmat(3,:);
    Rstruct(i,:) = R;
    %R = zeros(1,endcount)

    %shadows(i) = (1-sum(c))/sum(b);
    %abundance(i,:) = shadow*b+c;
    %reflectance(i,:) = c./abundance;
    %albedos = get_albedos(R,endmembers);
    %reconstruct = get_spectrum(b,c,albedos);

    %error = guassian_error_lim(inputmat,R,spectrum,endmembers);
end
R = mean(Rstruct);
albedos = get_albedos(R,endmembers);

for i=1:pixelcount
    spectrum = pixels(i,:);
    inputmat = fmincon(@(I) guassian_error_lim(I,R,spectrum,endmembers), ...
                default(1:2,:),[],[],[],[],zeros(2,endcount),ones(2,endcount),[], ...
                optimset('MaxFunEvals',6000, 'Display','off', ...
                'algorithm','interior-point'));
    b = inputmat(1,:);
    c = inputmat(2,:);
    %R = zeros(1,endcount)

    shadows(i) = (1-sum(c))/sum(b);
    abundance(i,:) = shadows(i)*b+c;
    reflectance(i,:) = c./abundance(i,:);
    P(i) = sum(c);
    reconstruct(i,:) = get_spectrum(b,c,albedos);

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


function [error] = guassian_error_lim(inputmat, R, spectrum, endmembers)
b = inputmat(1,:);
c = inputmat(2,:);

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
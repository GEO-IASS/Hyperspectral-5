function [ abundance, P, reflectance, reconstruct, error ] = multilin_custom(pixel, endmembers,options)

[endcount,~] = size(endmembers);
default = zeros(2,endcount);
default(2,:) = SCLSU(endmembers,pixel);

minimum = zeros(2,endcount);

if options.freeP
    
    minimum(1,:) = -1;
end

[inputmat,error] = fmincon(@(I) calc_error(I,pixel,endmembers,options), ...
            default(:,:),[],[],[],[],minimum,ones(2,endcount),[], ...
            optimset('MaxFunEvals',options.accuracy, 'Display','off', ...
            'algorithm','interior-point'));

[reflectance, abundance] = calc_param(inputmat,options);
reconstruct = build_spectrum(reflectance, abundance,endmembers);
P = sum(abundance.*reflectance);



end


function [P,a] = calc_param(inputmat,options)
P = inputmat(1,:);
a = inputmat(2,:);

if true %options.sumone
    a = a/sum(a);
end

if options.indepp
    P = P*0+mean(P);
end


end

function error = calc_error(inputmat,pixel,endmembers,options)

[P,a] = calc_param(inputmat,options);
error = get_error(pixel,P,a,endmembers);


end

function spec = build_spectrum(P,a,endmembers)
spec = (((1-P).*a)*endmembers)./(1-(P.*a)*endmembers);
end

function error = get_error(pixel,P,a,endmembers)
error = sumsqr(pixel-build_spectrum(P,a,endmembers));
end
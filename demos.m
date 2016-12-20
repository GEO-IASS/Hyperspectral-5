function s = demos
s.mlin = @multilinetest;
s.lib = @libtest;
s.alina = @Alinatest;
s.ac = @Alina_compare;
end


function [linerr,semerr] = Alina_compare()
close all;
tic()
[libs,data,~,rl] = load_Alina();
[~,semerr,linerr,P] = library_fast(libs,rl,data);
disp('making figures')
%close all;
figure;

density = reshape(log(linerr./semerr),13,19);
image(density,'CDataMapping','scaled');
colorbar;
colormap('jet')
print('Relative_Error','-dpdf')

figure;

density = reshape(log(semerr./linerr),13,19);
image(density,'CDataMapping','scaled');
colorbar;
colormap('jet')
print('Relative_Error_i','-dpdf')


figure;
density = reshape(P,13,19);
image(density,'CDataMapping','scaled');
colorbar;
colormap('jet')
caxis([0 1])
print('P','-dpdf')



disp(['Time elapsed: ' datestr(datenum(0,0,0,0,0,toc()),'HH:MM:SS')])
end

function Alinatest(lin)
close all;
tic()
[libs,data] = load_Alina();
[~, abundance, P, reflectance, ...
    reconstruct, albedos, error] = libraryperm(data,libs,true,lin,false,false);
disp('making figures')
%close all;
figure;
for i=1:4;
    subplot(2,2,i);
    density = reshape(abundance(:,i),13,19);
    image(density,'CDataMapping','scaled');
    colorbar;
    colormap('jet')
    caxis([0 1])
end
if lin
    print('linMESMA','-dpdf')
else
    print('semilinMESMA','-dpdf')
end
% figure;
% density = reshape(P(:),13,19);
% image(density,'CDataMapping','scaled');
% colorbar;
% colormap('jet')
% caxis([0 1])
% figure;
% density = reshape(sum(abundance(:,:).'> 0),13,19);
% image(density,'CDataMapping','scaled');
% colorbar;
% colormap('jet')
% 
% 



disp(['Time elapsed: ' datestr(datenum(0,0,0,0,0,toc()),'HH:MM:SS')])
end

function multilinetest(pixs,bands,ends)

[pixels, endmem, sol] = build_test_spectrum(pixs,bands,ends);
disp(sol)

[ abundance, reflectance, reconstruct,albedos, error ] = multilin(pixels,endmem);
disp(error)

%[specsort, X] = sort(spec);
close all;
figure;
hold on;
for i=1:pixs;
    plot(reconstruct(i,:),'b-p');
    plot(pixels(i,:),'r-h');
end


figure;
hold on;

boxplot((reconstruct-pixels).');

figure;
hold on; 
for i=1:ends
   plot(endmem(i,:),'k:x')
   %plot(albedos(i,:),'g:')
    
end

end

function libtest(libs,pixs,bands,ends,libsize,noise)


corr = randi([1 libsize],1,ends);
[pixels, library, sol] = build_test_lib(pixs,bands,libs,libsize,corr,noise);
disp(sol)

[sub_indices, abundance, P, reflectance, ...
    reconstruct, albedos, error] = librarymatch(pixels,library,1,Inf());
disp(error)
disp([num2str(corr) ' == ' num2str(sub_indices)])

%[specsort, X] = sort(spec);
close all;
figure;
hold on;

for i=1:pixs;
    plot(reconstruct(i,:),'b-p');
    plot(pixels(i,:),'r-h');
    plot(sol.pix(i,:),'g-s');
end
xlabel('band')
ylabel('intensiteit')
title('spectra van elke pixel')
legend('reconstructie','gegeven spectrum','exact spectrum')

figure;
hold on;
for i=1:pixs;
    plot(reconstruct(i,:) - pixels(i,:),'r-h');
    plot(reconstruct(i,:) - sol.pix(i,:),'g-s');
end

title('reconstructieerror op spectra van elke pixel')
legend('tov gegeven spectrum','tov exact spectrum')

xlabel('band')
ylabel('error')
figure;

boxplot((reconstruct-pixels).');
title('reconstructieerror tov gegeven spectrum')
xlabel('pixel')
ylabel('error')

figure;

boxplot((reconstruct-sol.pix).');
title('reconstructieerror tov exact spectrum')
xlabel('pixel')
ylabel('error')

figure;
hold on;
title('error op albedos')
for i=1:size(albedos,1);
    plot(albedos(i,:)-sol.w(i,:),'r-+');
end
xlabel('band')
ylabel('verschil albedos')

figure;

hold on;
subplot(2,2,1)
bar(abundance,'stacked')
title('abundanties')
xlabel('pixel')
ylabel('abundantie')
axis([0.5 pixs+0.5 0 1])

subplot(2,2,2)
hold on;
bar(sol.a,'stacked')
title('exacte abundanties')
xlabel('pixel')
ylabel('abundantie')
axis([0.5 pixs+0.5 0 1])

subplot(2,2,3)
bar(reflectance,'stacked')
title('reflectiekans')
xlabel('pixel')
ylabel('reflectiekans')
axis([0.5 pixs+0.5 0 Inf])

subplot(2,2,4)
hold on;
bar(sol.P,'stacked')
title('exacte reflectiekans')
xlabel('pixel')
ylabel('reflectiekans')
axis([0.5 pixs+0.5 0 Inf])

end
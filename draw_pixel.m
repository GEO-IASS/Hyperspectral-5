function d = draw_pixel(dataset,x,y)
waves = 0.4:0.0094594594594595:2.5;
waves([104:109,150:164,220]) = [];
spectrum = permute(dataset(x,y,:),[3 1 2]);
waves(end)
length(waves)
length(spectrum)

plot(waves,spectrum)
title('Spectrum')
xlabel('golflengte (nm)')
ylabel('Intensiteit')

 

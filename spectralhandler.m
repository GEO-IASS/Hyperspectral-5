function s = spectralhandler
s.draw_pixel = @draw_pixel;
s.display = @display_band;
s.slider = @slider;
s.PPIU = @simple_PPIU;
s.PPILH = @simple_PPILH;
s.sim = @sim;
s.show = @slider_density;
s.show_dataset = @show_dataset;
s.load_indian_pines = @load_indian_pines;
s.DEBUG_finish = @finish;
end

%% Datasets

function data = load_indian_pines()
load Indian_pines_corrected
data.intensity = indian_pines_corrected;
data.waves = 400:9.4594594594595:2500;
data.waves([104:109,150:164,220]) = [];
end

%% Display

function draw_pixel(dataset,x,y)
spectrum = get_pixel_spectrum(dataset.intensity,x,y);
draw_spectrum(spectrum,dataset.waves)
end

function p = draw_spectrum(spectrum,waves)
p = plot(waves,spectrum,'k');
title('Spectrum')
xlabel('golflengte (nm)')
ylabel('Intensiteit')

xlim([min(waves) max(waves)])
ylim([0 10000])
patch([620 620 750 750],[1000 2000 2000 1000],'r')
patch([495 495 620 620],[1000 2000 2000 1000],'g')
patch([450 450 495 495],[1000 2000 2000 1000],'b')
set(gca,'children',flipud(get(gca,'children')))
end

function p = draw_spectrum_custom(spectrum,waves,color)
p = plot(waves,spectrum,color);
title('Spectrum')
xlabel('golflengte (nm)')
ylabel('Intensiteit')

xlim([min(waves) max(waves)])
ylim([0 10000])
patch([620 620 750 750],[1000 2000 2000 1000],'r')
patch([495 495 620 620],[1000 2000 2000 1000],'g')
patch([450 450 495 495],[1000 2000 2000 1000],'b')
set(gca,'children',flipud(get(gca,'children')))
end

function plot = show_dataset(dataset)
red = get_band(dataset,560);
green = get_band(dataset,510);
blue = get_band(dataset,475);
plot = image(cat(3, stretch(red), stretch(green), stretch(blue)));
axis square;
end

function slider(dataset)
x = 1:10;
hplot = image(get_band(dataset,400),'CDataMapping','scaled');
h = uicontrol('style','slider','units','pixel','position',[20 0 200 20]);
addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,x,hplot,dataset));
end

function makeplot(hObject,event,x,hplot,dataset)
wavelen = get(hObject,'Value')*(2500-400)+400;
set(hplot,'CData',get_band(dataset,wavelen));
title(['Intensiteit bij ' num2str(correct_wavelength(wavelen)) ' nm'])
drawnow;
end

function slider_density(solution)

[endm bands] = size(solution.endmembers);
index = 1;
subplot(1,2,1)
plot_abundance = image(get_density_band(solution.density,index),'CDataMapping','scaled');
colorbar
colormap('hot')
%caxis([0 1])
axis square;
title(['band ' int2str(index)])
subplot(1,2,2)
plot_spectrum = draw_spectrum(solution.endmembers(index,:),solution.waves);
hold on;
title('spectrum')
h = uicontrol('style','slider','position',[20 0 500 20],'Min',0,'Max',endm,'SliderStep',[1/(endm+1) 0.5],'Value',2);
addlistener(h,'ContinuousValueChange',@(hObject, event) makeplot_density(hObject, event,plot_abundance,plot_spectrum,solution));
hold off;

end

function makeplot_density(hObject,event,plot_abundance,plot_spectrum,solution)

index = ceil(get(hObject,'Value')-0.5);
if index > 0
    set(plot_abundance,'CData',get_density_band(solution.density,index));
    set(plot_spectrum,'Ydata',solution.endmembers(index,:));
    subplot(1,2,1)
    title(['spectrum ' num2str(index)])
else
    set(plot_abundance,'CData',solution.error);
    set(plot_spectrum,'Ydata',solution.endmembers(1,:)*0);
    subplot(1,2,1)
    title('error ')
end
drawnow;
hold off;
end



%% tensor extraction and wavelenght selection
function data = get_band(dataset,wavelength)
bandindex = get_bandindex(dataset,wavelength);
data = dataset.intensity(:,:,bandindex);
end

function data = get_density_band(dataset,index)
data = dataset(:,:,index);
end

function spec = get_pixel_spectrum(dataset,x,y)
spec = permute(dataset(x,y,:),[3 1 2]);
end

function display_band(dataset,wavelength)
realwavelen = correct_wavelength(wavelength);
data = get_band(dataset,wavelength);
image(data,'CDataMapping','scaled');
colorbar
axis square 
title(['Intensiteit bij ' num2str(realwavelen) ' nm'])
xlabel('x')
ylabel('y')
end

function waves = waves(mini,maxi,steps,missing)
waves = 400:9.4594594594595:2500;
waves([104:109,150:164,220]) = [];
end

function bandindex = get_bandindex(dataset,wavelength)
diff = abs(dataset.waves - wavelength);
[c idx] = min(diff);
bandindex = idx;
end

function c = correct_wavelength(w,dataset)
index = get_bandindex(w);
wav = dataset.waves;
c = wav(index);
end

%% Vector and matrix handling

function s=sim(vector1,vector2)
s= norm(vector1/norm(vector1) - vector2/norm(vector2));
end

function s = stretch(matr)
s = (matr - min(min(matr)))/(max(max(matr))-min(min(matr)));
end

%% PPI

function p = PPI(dataset,number_of_endmembers,mean_hits)
%Runs a combination of the Pixel purity algoritm and the unconstained
%inversion algoritm
% this returns a strcuture with following properties:
% .dataset -> the original dataset w*h*b
% .endmembers -> a matrix containing the different endmembers as rows b*e
% .errors -> a matrix containing the error w*h
% .density -> a tensor containing the abundances of the endmembers w*h*e
[w h bands] = size(dataset.intensity);
Purity = sparse(w,h);
number_of_iterations = floor(w*h*mean_hits);

for i= 1:number_of_iterations
    vector = rand(bands,1)-0.5; % create random vector
    vector = vector / norm(vector); % normalise vector
    local = zeros([w,h]); % initalise local matrix
    %local contains the projection lenghts of each pixel's spectrum to the
    %above vector
    I = dataset.intensity;
    parfor x=1:w
        for y=1:h
            spec = get_pixel_spectrum(I,x,y); %get spectrum of pixel x,y
            local(x,y) = dot(spec,vector);%/norm(spec); %set projection lenght
        end
    end
    
    [num] = max(local(:)); % find maximum projection lenght
	[x y] = ind2sub(size(local),find(local==num)); % find indices of maximum projection lenght
    tx = x(1); %select one point 
    ty = y(1);
    [num] = min(local(:)); %alaog to above, but minimum this time
	[x y] = ind2sub(size(local),find(local==num));
    bx = x(1);
    by = y(1);
    
    Purity(tx,ty) = Purity(tx,ty)+1; %Add one to maximum
    Purity(bx,by) = Purity(bx,by)+1; %Add one to minimum
    if mod(i,10) == 0
        % This plots the whole basar
        subplot(1,2,2);
        image(local,'CDataMapping','scaled');
        title(['Index =' num2str(i) '/' num2str(number_of_iterations)]);
        colormap('hot')
        colorbar
        axis square;
        hold on;
        plot([ty by],[tx bx],'+k','linewidth',8)
        hold off;
        subplot(1,2,1);
        image(-Purity,'CDataMapping','scaled');
        caxis([-5 0])
        axis square;
        drawnow;
    end
end
endmembers = zeros(number_of_endmembers,bands);
%init endmember matrix
i = 1;
c = 1;
while i <= number_of_endmembers
    [num] = max(Purity(:)); %find pixel with max purity
	[x y] = ind2sub(size(Purity),find(Purity==num));
    x = x(1);
    y = y(1);
    spectrum = get_pixel_spectrum(dataset.intensity,x,y);%get spectrum of this pixel
    Purity(x,y) = 0;%removes this pixel so the next loop doesnt just reselect this
    %spectrum = spectrum / norm(spectrum);
    subplot(1,1,1)
    draw_spectrum(spectrum,dataset.waves)
    title(['Endmember ' num2str(i) '/' num2str(c)])
    drawnow;
    valid = 1;
    %for j=1:i-1
    %    if sim(endmembers(j,:)',spectrum) < 0.1 %Check if no other spectrum is very similar
    %        valid = 0;
    %        break
    %    end
    %end
    if valid == 1
        endmembers(i,:) = spectrum;%Add spectrum to endmember matrix
        i = i + 1;
    end
    c = c + 1;
    if c > 5000
        break
    end
end
p.endmembers = endmembers;% add endmembers to return object
p.intensity = dataset.intensity;
p.waves = dataset.waves;
end


%% Full Unmixing
function p = simple_PPIU(dataset,number_of_endmembers,mean_hits)
p = PPI(dataset,number_of_endmembers,mean_hits);
p.density = extract_unconstraned(p,number_of_endmembers);
finish(p);
end

function res = simple_PPILH(dataset,number_of_endmembers,mean_hits)
p = PPI(dataset,number_of_endmembers,mean_hits);
p.density = extract_Lawson_Hanson(p,number_of_endmembers);
res = finish(p);
end

function r = finish(dataholder)
[w h bands] = size(dataholder.density);
error = zeros(w,h);
for x=1:w %loop over pixels
    for y=1:h
        spec = get_pixel_spectrum(dataholder.intensity,x,y);
        alpha = permute(dataholder.density(x,y,:),[2 3 1]);
        recspec = dataholder.endmembers' * alpha';
        error(x,y) = norm(spec - recspec,2);
    end
end
dataholder.error = error;
slider_density(dataholder)
r = dataholder;
end

%% extraction

function density = extract_unconstraned(dataholder,number_of_endmembers)

[w h bands] = size(dataholder.intensity);

density = zeros(w,h,number_of_endmembers);

extract = pinv(dataholder.endmembers');
for x=1:w %loop over pixels
    parfor y=1:h
        spec = get_pixel_spectrum(dataholder.intensity,x,y);% Get spectrum of pixel x,y
        alpha = extract * spec;
        % multiply spectrum with Moore-Penrose inverse
        alpha = permute(alpha,[1 3 2]);
        density(x,y,:) = alpha;
    end
end

end

function density = extract_Lawson_Hanson(dataholder,number_of_endmembers)

[w h bands] = size(dataholder.intensity);

density = zeros(w,h,number_of_endmembers);

for x=1:w %loop over pixels
    for y=1:h
        spec = get_pixel_spectrum(dataholder.intensity,x,y);% Get spectrum of pixel x,y
        alpha = lsqnonneg(dataholder.endmembers',spec);
        alpha = permute(alpha,[1 3 2]);
        density(x,y,:) = alpha;
        draw_spectrum_custom(spec,dataholder.waves,'r')
        
        hold on;
        specr = spec*0;
        for e=1:number_of_endmembers
            espec = alpha(e) * dataholder.endmembers(e,:);
            specr = specr + espec';
            draw_spectrum(espec,dataholder.waves)
            
        end
        
        draw_spectrum_custom(specr,dataholder.waves,'b');
        title(['position (' num2str(x) ',' num2str(y) ')']);
        hold off;
        drawnow;
        
    end
end

end





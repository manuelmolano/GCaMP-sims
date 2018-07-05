function [final_neuron_shape,index_area_comp_unif,actual_sigma_tam_neuronas] = build_neuron(binornd_neuronas,tam,sigma_tam_neuronas,nucleus_prop,posicion_neuronas,...
    nucleus_weight,kernel_smoothing_neurons,neuron_threshold,neuron_homogeneity)
[X,Y] = meshgrid(1:tam);

%sigmas for gaussian that define the donut shape of the neuron
actual_sigma_tam_neuronas = sigma_tam_neuronas + randn*sigma_tam_neuronas/20;
sigma_tam_nucleus = actual_sigma_tam_neuronas*nucleus_prop;
%gaussian distribution
positive_gaussian = (1/(2*pi*actual_sigma_tam_neuronas^2))*...
    exp(-((X-posicion_neuronas(1)).^2)/(2*actual_sigma_tam_neuronas^2)-((Y-posicion_neuronas(2)).^2)/(2*actual_sigma_tam_neuronas^2));
negative_gaussian = (1/(2*pi*sigma_tam_nucleus^2))*...
    exp(-((X-posicion_neuronas(1)).^2)/(2*sigma_tam_nucleus^2)-((Y-posicion_neuronas(2)).^2)/(2*sigma_tam_nucleus^2));

positive_gaussian = positive_gaussian/sum(positive_gaussian(:));
negative_gaussian = negative_gaussian/sum(negative_gaussian(:));
%subtract the nucleus
rng(binornd_neuronas);
original_shape = positive_gaussian - min(1,max(0,randn*0.25+nucleus_weight))*negative_gaussian;
original_shape(original_shape>1) = 1;
original_shape(original_shape<0) = 0;
%get the irregular shape (irregularity is controled by neuron_homogeneity
aux = neuron_homogeneity*original_shape/max(original_shape(:));
aux(aux>1) = 1;
rng(binornd_neuronas);
neuron_shape = binornd(1,aux);%this is my donut-shaped neuron
%median filter
neuron_shape_filtered = medfilt2d(neuron_shape,1);
%control if the median filter has removed all pixels
counter_aux = 0;
while nnz(neuron_shape_filtered)==0
    neuron_shape = binornd(1+counter_aux,aux);
    neuron_shape_filtered = medfilt2d(neuron_shape,1);
    counter_aux = counter_aux + 1;
end
%convolutions
neuron_shape_smoothed = conv2(neuron_shape_filtered,kernel_smoothing_neurons,'same');
final_neuron_shape = neuron_shape_smoothed;%.*neuron_shape_filtered;
final_neuron_shape = final_neuron_shape/max(final_neuron_shape(:));
final_neuron_shape(final_neuron_shape<neuron_threshold) = 0;
max_distance = 18.5;%pixels
[X1, X2] = meshgrid(1:tam);
pixels_patch = [X1(:) X2(:)];
index_area_comp_unif = (pixels_patch(:,1)-posicion_neuronas(1)).^2 + (pixels_patch(:,2)-posicion_neuronas(2)).^2>(max_distance/2)^2;
final_neuron_shape(index_area_comp_unif) = 0;
index_area_comp_unif = find(index_area_comp_unif==0);

%{
    figure
    subplot(2,3,1)
    imagesc(original_shape)
    axis image
    title('original shape')
    subplot(2,3,2)
    imagesc(neuron_shape)
    axis image
    title('effective shape')
    subplot(2,3,3)
    imagesc(neuron_shape_filtered)
    axis image
    title('filtered shape')
    subplot(2,3,4)
    imagesc(neuron_shape_smoothed)
    axis image
    title('smoothed')
    subplot(2,3,5)
    imagesc(final_neuron_shape)
    axis image
    title(['final shape. Num pixels=' num2str(nnz(final_neuron_shape))])
    
   %} 

end


function filtered_image = medfilt2d(image,size_w)
%if a pixel is completely isolated within a vicinity defined by size_w it is made 0
[i, j] = find(image);

%crop the image to avoid going over all pixels (most of which will be 0)
crop_image = image(min(i)-size_w:max(i)+size_w,min(j)-size_w:max(j)+size_w);
crop_filtered_image = zeros(size(crop_image));
for ind1=size_w+1:size(crop_image,1)-size_w
    for ind2=size_w+1:size(crop_image,2)-size_w
        window = crop_image(ind1-size_w:ind1+size_w,ind2-size_w:ind2+size_w);
        if window(size_w+1,size_w+1)~=0 && nnz(window)>1
            crop_filtered_image(ind1,ind2) = 1;
        end
    end
end

filtered_image = zeros(size(image));

filtered_image(min(i)-size_w:max(i)+size_w,min(j)-size_w:max(j)+size_w) = crop_filtered_image;
end

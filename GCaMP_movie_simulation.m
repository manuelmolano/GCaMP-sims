
function data = ...
    GCaMP_movie_simulation(tam,duration_ms,dwellTime,num_neurons,neuron_homogeneity,neuron_threshold,nucleus_weight,sigma_tam_neuronas,...
    spike_rate,baseline_rate,noise_photon_level,neuropil_level,rise_time,decay_time_constant,pnonlin,transient_peak,...
    mat_surrounding,num_pixels_th,folder,name, pintar, pintar_all)
%this function creates a toy movie in which neurons get activated acording
%to an inegrate and fire model adapted from the material of the Theoretical Neuroscience course (coursera).
%It takes several input parameters:
%data: a struct where to save the movie
%tam (size of the image)
%sigma_tam_neuronas (controls neurons size)
%spike_rate (spike rate of all neurons)
%num_neurons, which should be an even number because neurons' activity is created in pairs (to simulate synchrony)
%duration_ms, duration of the movie
%syncrho (neurons are sinchronize on a pair-wise fashion following this paramater; 0 --> no synchro, 1 --> perfect synchro)
%dwellTime (s), time the laser spends on each pixel
%pintar: whether to plot or not the figures
%noise_photon_level (controls the variance of the uniform noise we are adding to the movie)
close all
global sampling_noise_constant
sampling_noise_constant = 0.1;
visible = 'on';
colores = [[1 0 0];[0 1 0];[0 0 1];[0 1 1];[1 0 1];[1 1 0];[0 0 0]];%colores is the variable that keeps the colours used to draw the ROIs
colores = [colores(2:end,:);colores(1:end-1,:)*0.9;colores(1:end-1,:)*0.8,;colores(1:end-1,:)*0.7...
    ;colores(1:end-1,:)*0.6;colores(1:end-1,:)*0.5,;colores(1:end-1,:)*0.4...
    ;colores(1:end-1,:)*0.3;colores(1:end-1,:)*0.2,;colores(1:end-1,:)*0.1];%there are up to 25 different colours
colores = [colores;colores];
mat_surrounding = [-1,mat_surrounding];
%saving info
disp(name)
name(strfind(name,'.')) = '';
folder_name = [folder name '/'];
if ~exist(folder_name,'dir')
    mkdir(folder_name)
end


%value for isolated spikes
%these are the minimum interspike distances (pre and post)
baseline_window = 1000;
peak_window = 400;
%stim times. If greater than 0 stim_times will be taken as the stimulus onset
baseline_used_for_exps = 7000;
response_used_for_exps = 2000;
baseline_duration = 0.25;
response_duration = baseline_duration/4;
stim_times = round((baseline_duration:response_duration+baseline_duration:1-response_duration)*duration_ms);
baseline_duration = round(baseline_duration*duration_ms);
response_duration = (response_duration*duration_ms);
assert(baseline_duration>baseline_used_for_exps)
assert(response_duration>response_used_for_exps)
%neurons' position
posicion_neuronas =...
    mosaico_random_induccion(tam,sigma_tam_neuronas,num_neurons);%this function gives the positions of the neurons, following some contraints.
num_neurons = size(posicion_neuronas,2);

%get spike trains.
%here I obtain the spike train of each neuron
spikes = zeros(num_neurons,duration_ms);
fr_resp = zeros(num_neurons,50);
fr_bl = zeros(num_neurons,1);
for ind_neuronas=1:num_neurons
    [spikes(ind_neuronas,:),fr_resp(ind_neuronas,:),fr_bl(ind_neuronas)] =...
        circuit_II(duration_ms, spike_rate, stim_times, baseline_rate);
end

%get the frame period and frame times of the raster-scan movie
frame_period_ground_truth = 0.001;%s this is to simulate the image adquisition rate
pixels_per_line = tam;
linesPerFrame = tam;
framePeriod = dwellTime*pixels_per_line*linesPerFrame;
frame_times =  0:framePeriod:duration_ms/1000-framePeriod;
average_window_choose_roi = 0.4;%this is the period (s) that will be taken to average frames around an event
average_window_choose_roi_current = ceil(average_window_choose_roi/framePeriod);%here we use the frame period to infer the number of frames
neuropil_homogeneity = 0.9;%controls the uniformity of the neuropil (1 --> uniform)
neuropil_baseline_level = 0.5;
nucleus_prop = 0.9;%defines the relative size of the nucleus
neurons_baseline_level = 0.1;%%contribution of neuron's baseline
%smoothing kernel
index = -1:0.2:1;
smoothing_factor = 0.05;%we will convolve the frames with a square pulse to remove some salt and pepper noise. This is the size of the pulse.
kernel_smoothing_neurons = smoothing_kernel(index, smoothing_factor);
%smoothing kernel
smoothing_factor = 1;%we will convolve the frames with a square pulse to remove some salt and pepper noise. This is the size of the pulse.
t_trace = 1:duration_ms;%this is the fine-resolution time that I use for the neurons activity

%get rid of non-active neurons. The number of neurons is updated. [having
%non-active cells is an issue because the function get_rois_and_activities does not find
%pixels above baseline in line: [i,j] = find(reference_aux>0);]
%first remove spikes that are too close to the end
spikes(:,end-ceil(2*framePeriod*1000):end) = 0;
total_activity = sum(spikes,2);
spikes(total_activity==0,:) = [];
posicion_neuronas(:,total_activity==0) = [];
num_neurons = size(spikes,1);

%here I am creating the random seed for each neuron, so their shapes is the
%same everytime I add them to the image (when they fire)
binornd_neuronas = cell(1,num_neurons + 1);
for ind_n=1:num_neurons+1
    binornd_neuronas{ind_n} = rng(ind_n);
end

%obtain the baseline contribution for each neuron (how much fluorescent when not firing)
baseline_contr = rand(1,num_neurons);



%this matrix keeps the activity of each neurons along the duration of the movie
mat_fluorescence = zeros(num_neurons+1,duration_ms);
%these are samples for the imaging (depend on the ground truth frame-rate)
samples = round(1:round(1000*frame_period_ground_truth):duration_ms);
movie = zeros(numel(samples),tam^2);%the movie without noise
movie_baseline = zeros(1,tam^2);%the movie without noise
movie_neuropil = zeros(numel(samples),tam^2);%the movie without noise
ground_truth_spiketimes = cell(1,num_neurons+1);
pixels_per_neuron =  cell(1,num_neurons);
actual_sigma_tam_neuronas = zeros(1,num_neurons);
total_num_pixels_per_neuron = zeros(1,num_neurons);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get neuropil
[r,patch_pixels_original] = neuropil_activity(duration_ms);
for ind_n=1:num_neurons
    patch_pixels = patch_pixels_original+repmat(posicion_neuronas(end:-1:1,ind_n)',size(patch_pixels_original,1),1);
    index = sub2ind([tam tam],patch_pixels(:,1),patch_pixels(:,2));
    aux = r/std(r(:));
    movie_neuropil(:,index) = aux;
end

movie_neuropil = movie_neuropil*neuropil_level;
movie_neuropil = movie_neuropil-min(movie_neuropil(:));

%for each neuron, I create its fluorescence-activity profile based on the spikes I generated before
for ind_n=1:num_neurons
    %get spikes
    spike_times = find(spikes(ind_n,:));
    ground_truth_spiketimes{ind_n} = spike_times;
    %Here I put the transients in different rows of g_syn and the sum them
    g_syn = zeros(numel(spike_times),duration_ms);%alpha func synaptic conductance
    for ind_spk=1:numel(spike_times)
        %get decay trace
        trace = max(0,min(1,(t_trace-spike_times(ind_spk)-rise_time))).*...
            exp(-max(0,(t_trace-spike_times(ind_spk)-rise_time))/decay_time_constant);
        %get rise period
        rise_period = spike_times(ind_spk):spike_times(ind_spk)+rise_time;
        rise_phase = 0:numel(rise_period)-1;
        rise_phase = rise_phase/max(rise_phase);
        trace(rise_period) = rise_phase;
        %apply non-linearity as in Akerboom et al.
        trace = polyval([fliplr(pnonlin) 1 0],trace);
        g_syn(ind_spk,:) = transient_peak*trace;
        %{
        transient = g_syn(ind_spk,:);
        %plot calcium transient
        tam_screen = get(0,'ScreenSize');
        index = find(transient>=max(transient)-max(transient)/20);
        numel(index)
        figure('Position',[10 10 tam_screen(3)-50 tam_screen(4)-50]);
        plot(transient,'lineWidth',3)
        hold on
        plot(index,transient(index),'lineWidth',3)
        ylim([-0.01 0.11])
        numel(index)
        xlabel('time (ms)')
        ylabel('fluorescence (a.u.)')
        set(gca,'fontSize',20)
        legend('transient','within the 5% of the peak')
        axes_rasters = axes('Position',[.5 .5 .3 .3]);
        plot(axes_rasters,1810:1870,transient(1810:1870),'lineWidth',3)
        hold on
        plot(axes_rasters,index,transient(index),'lineWidth',3)
        ylim([0.095 0.101])
        xlabel('time (ms)')
        ylabel('fluorescence (a.u.)')
        set(gca,'fontSize',16)
        %}
    end
    mat_fluorescence(ind_n,:) = sum(g_syn,1);%this is the fluorescence-activity of neuron ind_n
    assert(nnz(mat_fluorescence(ind_n,:)~=0)>0)
    

    %get neuron shape
    [final_neuron_shape,~,actual_sigma_tam_neuronas(ind_n)] = ...
        build_neuron(binornd_neuronas{ind_n},tam,sigma_tam_neuronas,nucleus_prop,posicion_neuronas(:,ind_n),...
    nucleus_weight,kernel_smoothing_neurons,neuron_threshold,neuron_homogeneity);
    
    
    index_n = find(final_neuron_shape);
    total_num_pixels_per_neuron(ind_n) = numel(index_n);
    [i,j] = find(final_neuron_shape);
    pixels_per_neuron{ind_n} = [j,i];
   
    
    %once I get the neuron shape, at each frame I multiply it by the fluorescence
    for ind_frame=samples
        movie(ind_frame,:) = movie(ind_frame,:) + final_neuron_shape(:)'*mat_fluorescence(ind_n,ind_frame);
    end
    movie_baseline = movie_baseline + baseline_contr(ind_n)*final_neuron_shape(:)';
    
end


%here I am adding the mean to have basal flurescence on each frame
movie = movie + neurons_baseline_level*repmat(movie_baseline,size(movie,1),1);
movie_final = movie +  movie_neuropil;


if pintar
    tam_screen = get(0,'ScreenSize');
    %main summary figure
    path_panel_num_pixels = mat_surrounding(max(1,round(numel(mat_surrounding)/4)));
    %summary figure
    h_paper = figure('name','paper','OuterPosition',[10 10 tam_screen(3) tam_screen(4)]);
    [axes_rasters,axes_stats] = figure_panels();
    axes(axes_rasters(1))
    frame = randi(size(movie,1)-1000)+1000;
    aux = movie(frame,:);
    aux = reshape(aux,tam,[]);
    aux(size(aux,1)-2:size(aux,1),1:50) = max(aux(:));
    imagesc(aux)
    colormap('gray')
    axis image
    title('neural activity')
    set(axes_rasters(1),'fontSize',20)
    axes(axes_rasters(3))
    aux = movie_neuropil(frame,:);
    aux = reshape(aux,tam,[]);
    imagesc(aux)
    colormap('gray')
    axis image
    title('neuropil activity')
    set(axes_rasters(3),'fontSize',20)
end

if pintar_all
    figure('units','normalized','position',[0.1 0.1 0.5 0.5]);
    hist(total_num_pixels_per_neuron)
    xlabel('number of pixels')
    ylabel('number of neurons')
    set(gca,'fontSize',20)
    axis square
    title(['mean: ' num2str(mean(total_num_pixels_per_neuron)) '(' num2str(std(total_num_pixels_per_neuron)) ')'])
    %plot histogram of number of pixels
    
    figure('units','normalized','position',[0.1 0.1 0.5 0.5]);
    imagesc(mat_fluorescence)
    set(gca,'xTick',[],'yTick',[])
    xlabel('time')
    ylabel('neuron')
    set(gca,'fontSize',20)
    axis square
    
    %play movie
    h8 = figure('name','average','OuterPosition',[10 10 tam_screen(3)/2 tam_screen(4)/2]);
    for ind=1:100:size(movie,1)
        subplot(2,2,1)
        aux = movie(ind,:);
        aux = reshape(aux,tam,[]);
        imagesc(aux)%,[0 max(movie(:))])
        colormap('gray')
        axis image
        title(['neurons activity (frame=' num2str(ind) ')'])
        subplot(2,2,2)
        aux = movie_neuropil(ind,:);
        aux = reshape(aux,tam,[]);
        imagesc(aux)%,[0 max(movie_neuropil(:))])
        colormap('gray')
        axis image
        title('neuropil activity')
        subplot(2,2,4)
        aux = movie_final(ind,:);
        aux = reshape(aux,tam,[]);
        imagesc(aux)%,[0 max(movie_final(:))])
        colormap('gray')
        axis image
        title('high resolution activity (no sampling noise)')
        getframe;
    end
end

%put stuff on the movie_doc struct
data_toy.average_movie = mean(movie,1);
data_toy.binornd_neuronas = binornd_neuronas;
data_toy.actual_sigma_tam_neuronas = actual_sigma_tam_neuronas;
data_toy.position = posicion_neuronas;
data_toy.tam = tam;%microns
data_toy.linesPerFrame = tam;
data_toy.pixels_per_line = tam;
data_toy.sigma_tam_neuronas = sigma_tam_neuronas;%microns
data_toy.spike_rate = spike_rate;%Hz
data_toy.baseline_rate = baseline_rate;%Hz
data_toy.num_neurons = num_neurons;%num neurons
data_toy.duration_ms = duration_ms;%ms
data_toy.dwellTime = dwellTime;
data_toy.framePeriod = framePeriod;
data_toy.frame_times =  frame_times;
data_toy.average_window_choose_roi = average_window_choose_roi;%this is the period (s) that will be taken to average frames around an event
data_toy.average_window_choose_roi_current = average_window_choose_roi_current;%here we use the frame period to infer the number of frames
data_toy.neuropil_homogeneity = neuropil_homogeneity;%controls the uniformity of the neuropil (1 --> uniform)
data_toy.neuron_homogeneity = neuron_homogeneity;%controls the shape of neurons through the gaussian distribution (1 --> shape = gaussian distr.)
data_toy.neuropil_level = neuropil_level;
data_toy.noise_photon_level = noise_photon_level;
data_toy.neuropil_baseline_level = neuropil_baseline_level;
data_toy.neurons_baseline_level = neurons_baseline_level;
data_toy.smoothing_factor = smoothing_factor;%we will convolve the frames with a square pulse to remove some salt and pepper noise. This is the size of the pulse.
data_toy.frame_rate_ground_truth = frame_period_ground_truth;%this is to simulate the image adquisition rate
data_toy.transient_peak = transient_peak;%this is the transient maximum
data_toy.rise_time = rise_time;%time (ms) it takes to reack the peak
data_toy.decay_time_constant = decay_time_constant; %this controls the decay time constant
data_toy.mat_surrounding = mat_surrounding;
data_toy.nucleus_prop = nucleus_prop;
data_toy.electro_data = ground_truth_spiketimes;
data_toy.stim_times = stim_times;
data_toy.baseline_duration = baseline_used_for_exps;
data_toy.response_duration = response_used_for_exps;
data_toy.nucleus_weight = nucleus_weight;
data_toy.total_num_pixels_per_neuron = total_num_pixels_per_neuron;
data_toy.num_pixels_th = num_pixels_th;
data_toy.pixels_per_neuron = pixels_per_neuron;


%now we sample the movie
%note that now we get the sampling time when we actually sample. We have a
%sampling times vector for each pixel and we then average those times for
%the pixels corresponding to the given neuron. Before, we were averaging
%those pixels and calculating the sampling times for the hypothetical
%centre pixel (it could actually happen that the centre pixel was not a
%real pixel. The method makes more sense because the sampling times are
%exaclty those used to sample the original movie.
[data.movie_doc.movie_ruido,raster_sampling_time] = movie_sampling_raster(movie_final,frame_times,dwellTime,tam,noise_photon_level);
data_toy.movie_doc.movie_ruido = data.movie_doc.movie_ruido;

if pintar
    %plot frame from raster movie in summary figure
    axes(axes_rasters(2))
    frame = randi(size(data.movie_doc.movie_ruido,1)-10)+10;
    aux = data.movie_doc.movie_ruido(frame,:);
    aux = reshape(aux,tam,[]);
    imagesc(aux')
    colormap('gray')
    axis image
    title('example frame')
end


%performance measurements
%correlation
mean_performance_raster = zeros(num_neurons,numel(mat_surrounding));
mean_performance_scan = zeros(num_neurons,numel(mat_surrounding));
%using stimulus response (Chen et al.)
SNRChen_raster = zeros(num_neurons,numel(mat_surrounding));
SNRChen_raster_max = zeros(num_neurons,numel(mat_surrounding));
SNRChen_raster_std = zeros(num_neurons,numel(mat_surrounding));
SNRChen_scan = zeros(num_neurons,numel(mat_surrounding));
SNRChen_max = zeros(num_neurons,numel(mat_surrounding));
SNRChen_std = zeros(num_neurons,numel(mat_surrounding));
SNRChen_gt = zeros(num_neurons,numel(mat_surrounding));
%from transient associated to isolated spikes
SNRMarco_scan = nan(num_neurons,numel(mat_surrounding));
SNRMarco_max = nan(num_neurons,numel(mat_surrounding));
SNRMarco_std = nan(num_neurons,numel(mat_surrounding));
frame_rates = zeros(1,numel(mat_surrounding));
num_pixels_total = zeros(1,numel(mat_surrounding));
if pintar_all
    %play a video of the recorded movie
    set(axes_rasters(2),'fontSize',20)
    h8 = figure('units','normalized','position',[0.3 0.3 0.3 0.5]);
    for ind=1:1:size(data.movie_doc.movie_ruido,1)
        aux = data.movie_doc.movie_ruido(ind,:);
        aux = reshape(aux,tam,[]);
        imagesc(aux')
        colormap('gray')
        axis image
        title(['neurons activity (frame=' num2str(ind) ')'])
        pause(tam*tam*dwellTime)
    end
    %define figures
    h7 = figure('name','individual neurons performance','OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
    h12 = figure('name','SNRChen','OuterPosition',[10 10 tam_screen(3) tam_screen(4)/2],'visible',visible);
    h13 = figure('name','SNRChen raster scan','OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
    h14 = figure('name','SNRChen smart-line scan','OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
    
end
panel_counter = 1;
%START SCANNING
%get the original path (without surrounding)
counter = 0;
scan_path_original = zeros(2,500*num_neurons);
tags = zeros(1,500*num_neurons);
mat_num_pixels_per_neurons = zeros(1, num_neurons);
scan_path_original_Lillis = zeros(2,500*num_neurons);
tags_Lillis = zeros(1,500*num_neurons);
mat_num_pixels_per_neurons_Lillis = zeros(1, num_neurons);
counter_Lillis = 0;
SNRGUI_curves = zeros(num_neurons,150);
for ind_n=1:num_neurons
    %get pixels in the vicinity of the neuron
    all_pixels_on_square = pixels_per_neuron{ind_n};%select_pixels(posicion_neuronas(:,ind_n), tam, data_toy.average_movie,sigma_tam_neuronas,sigma_factor);
    all_pixels_on_square_1D = sub2ind([tam tam],all_pixels_on_square(:,1),all_pixels_on_square(:,2));
    %find the strongest event and get the average trace around it
    imagen = find_strongest_event(data.movie_doc.movie_ruido,tam,all_pixels_on_square_1D,average_window_choose_roi_current,numel(frame_times));
    
    %select pixels maximizing the SNR
    [scan_path_original, tags, num_pixels,SNRGUI_curves(ind_n,:)] =...
        get_SNRGUI_for_different_num_of_pixels(imagen,all_pixels_on_square, all_pixels_on_square_1D, data.movie_doc.movie_ruido, ind_n, scan_path_original, num_pixels_th, counter, tags, pintar_all);
    
    [scan_path_original_Lillis, tags_Lillis, num_pixels_Lillis] =...
        get_SNRGUI_for_different_scanLines(imagen, all_pixels_on_square, ind_n, scan_path_original_Lillis, counter_Lillis, tags_Lillis, tam);
    assert(num_pixels_Lillis>0)
    counter_Lillis = counter_Lillis + num_pixels_Lillis;
    mat_num_pixels_per_neurons_Lillis(ind_n) = num_pixels_Lillis;
    assert(num_pixels>0)
    counter = counter + num_pixels;
    mat_num_pixels_per_neurons(ind_n) = num_pixels;
end
data_toy.num_pixels_per_neurons = mat_num_pixels_per_neurons;
data_toy.SNRGUI_curves = SNRGUI_curves;
scan_path_original = scan_path_original(:,1:counter);
tags = tags(:,1:counter);
scan_path_original_Lillis = scan_path_original_Lillis(:,1:counter_Lillis);
tags_Lillis = tags_Lillis(:,1:counter_Lillis);
data_toy.num_pixels_per_neurons_Lillis = mat_num_pixels_per_neurons_Lillis;

%here I calculate distances from all points to the scan trajectory
[aux1, aux2] =  meshgrid(1:tam,1:tam);
all_pixels = [aux2(:) aux1(:)];
scan_distances_to_trajectory = zeros(1,size(all_pixels,1));
scan_neuron_tag_all_pixels = zeros(1,size(all_pixels,1));
for ind_sc=1:size(all_pixels,1)
    [scan_distances_to_trajectory(ind_sc),index] =...
        min(sqrt((all_pixels(ind_sc,1)-scan_path_original(1,:)).^2 + (all_pixels(ind_sc,2)-scan_path_original(2,:)).^2));
    %here we get as well the tag of each pixels to quickly assign tags to any trajectory
    scan_neuron_tag_all_pixels(ind_sc) = tags(index);
end

for ind_surr=1:numel(mat_surrounding)
    if mat_surrounding(ind_surr)>=0
        %get extra pixels that surround the selected ones
        scan_extra_pixels = all_pixels(scan_distances_to_trajectory<=mat_surrounding(ind_surr),:);
        aux = sub2ind([tam tam],scan_extra_pixels(:,1),scan_extra_pixels(:,2));
        traj_tags = scan_neuron_tag_all_pixels(aux);
        %calculate the shortest path (genetic algorithm)
        scan_path = group_pixels_and_calculate_shortest_path(scan_extra_pixels,traj_tags)';
    elseif mat_surrounding(ind_surr)==-1
        %calculate the shortest path (genetic algorithm)
        scan_path = group_pixels_and_calculate_shortest_path_lineScan(scan_path_original_Lillis',tags_Lillis)';
    end
    
    %summary figure
    if pintar && mat_surrounding(ind_surr)==path_panel_num_pixels
        aux = mean(movie,1);
        aux = reshape(aux,tam,[]);
        imagesc(axes_rasters(4),aux)
        hold(axes_rasters(4),'on')
        plot(axes_rasters(4),posicion_neuronas(1,:),posicion_neuronas(2,:),'*g','markerSize',12)
        
    end
    if pintar_all
        %plot average projection and scan paths
        h1 = figure('OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
        aux = mean(movie,1);
        aux = reshape(aux,tam,[]);
        imagesc(aux')
        colormap('gray')
        axis image
        hold on
        plot(posicion_neuronas(2,:),posicion_neuronas(1,:),'*g','markerSize',12)
        plot(scan_path(2,:),scan_path(1,:),'+r')
        plot(scan_path(2,:),scan_path(1,:),'c')
    end
    
    % All pixels in the path. Need to include the pixels inbetween the
    % neurons, which is what the prairie software does (it cannot jump)
    whole_path = include_all_pixels(scan_path,tam,pintar_all);
    num_pixels_total(ind_surr) = size(whole_path,2);
    frame_rates(ind_surr) = 1/(dwellTime*size(whole_path,2));
    
    
    
    %sample the movie
    %note that now we get the sampling time when we actually sample. We have a
    %sampling_times vector for each pixel and we then average those times for
    %the pixels corresponding to the given neuron. Before, we were averaging
    %those pixels and calculating the sampling times for the hypothetical
    %centre pixel (it could actually happen that the centre pixel was not a
    %real pixel. The method makes more sense because the sampling times are
    %exaclty those used to sample the original movie.
    [movie_sample,movie_sample_only_neurons,scan_sampling_times] =...
        scan_path_movie(whole_path,dwellTime,data_toy,noise_photon_level,movie_final,movie);
    
    if pintar_all
        %plot activities from SmartLine Scan with and woithout noise
        h3 = figure('OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
        subplot(2,1,1)
        imagesc(movie_sample)
        subplot(2,1,2)
        plot(mean(movie_sample,1))
        xlim([1 size(movie_sample,2)])
        h4 = figure('OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
        subplot(2,1,1)
        imagesc(movie_sample_only_neurons)
        subplot(2,1,2)
        plot(mean(movie_sample_only_neurons,1))
        xlim([1 size(movie_sample_only_neurons,2)])
    end
    
    %get rois (segmentation) and activities from SmartLine scan
    try
        [activities_scan_path, frame_times_scan_path,tags,~, effective_pixels_per_neuron, effective_num_pixels_per_neuron] =...
            get_rois_and_activities(movie_sample,pixels_per_neuron, num_neurons,...
            whole_path, scan_sampling_times); %#ok<ASGLU>
    catch
        keyboard
    end
    %this has to be here because it uses the updated tags from get_rois_and_activities
    if pintar_all
        %plot whole scan path
        set(0,'currentFigure',h1)
        for ind_n=0:num_neurons
            plot(whole_path(2,tags==ind_n),whole_path(1,tags==ind_n),'*','color',colores(ind_n+1,:))
        end
        plot(whole_path(2,:),whole_path(1,:),'color',[.7 .7 .7])
    end
    if pintar && mat_surrounding(ind_surr)==path_panel_num_pixels
        %plot whole scan path also in the summary figure
        for ind_n=0:num_neurons
            plot(axes_rasters(4),whole_path(1,tags==ind_n),whole_path(2,tags==ind_n),'*','color',colores(ind_n+1,:))
        end
        title(axes_rasters(4),'example smart-line scan path')
        set(axes_rasters(4),'fontSize',20)
    end
    
    
    %get raster activities from the same pixels
    activities_raster = zeros(size(data.movie_doc.movie_ruido,1),num_neurons);
    timing_raster = zeros(size(data.movie_doc.movie_ruido,1),num_neurons);
    for ind_n=1:num_neurons
        %get the activities obtained from the raster
        [activities_raster(:,ind_n),timing_raster(:,ind_n)] =...
            get_raster_activities(data.movie_doc.movie_ruido,whole_path,tags,ind_n,tam,raster_sampling_time);
    end
    
    %get ground truth activities from the same pixels. The activities
    %extracted here will not have any delay between pixels because they are
    %obtained from the original movie. Thus, while for the raster and
    %scanLine, different pixels of a neuron are sampled at different times
    %(corresponding to different frames of the original movie), these
    %activities are the original ones with the activities of all pixels
    %corresponding to a neuron sampled at the same exact time. Note that
    %this can cause differences in the average signals, which correspond to
    %the same ms (the same original movie frame) in the case of the
    %ground truth activities but to an average across activities that
    %happened at sligthly different times for the scanLine and, especially,
    %the raster scan.
    activities_ground_truth = zeros(size(movie,1),num_neurons);
    for ind_n=1:num_neurons
        %get the activities obtained from the raster
        [activities_ground_truth(:,ind_n)] =...
            get_raster_activities(movie,fliplr(whole_path')',tags,ind_n,tam);
    end
    
    activities_ground_truth_noise_and_neuropil = zeros(size(movie,1),num_neurons);
    for ind_n=1:num_neurons
        %get the activities obtained from the raster
        [activities_ground_truth_noise_and_neuropil(:,ind_n)] =...
            get_raster_activities(movie_final,fliplr(whole_path')',tags,ind_n,tam);
    end
    
    %interpolate activities
    activities_raster_int = interpolate_activity(activities_raster,timing_raster,data_toy.duration_ms);
    activities_scan_int = interpolate_activity(activities_scan_path,frame_times_scan_path,data_toy.duration_ms);
    
    %measure performance
    %correlation between ground truth and scanned activities
    if pintar_all
        h6 = figure('OuterPosition',[10 10 tam_screen(3) tam_screen(4)],'visible',visible);
    end
    scan_corr = zeros(1,num_neurons);
    raster_corr = zeros(1,num_neurons);
    for ind_n=1:num_neurons
        ground_truth = activities_ground_truth(:,ind_n);
        ground_truth_noise_and_neuropil = activities_ground_truth_noise_and_neuropil(:,ind_n);
        %corr might not be the best quality measurement because it doesn't
        %provide an estimate of how well the sampled traces matches the
        %amplitude of the transients (which is important for spike
        %inference). I will try to implement a spike inference method that
        %would definitelly give us a good quality measurement.
        scan_corr(ind_n) = corr(activities_scan_int(:,ind_n),ground_truth);
        raster_corr(ind_n) = corr(activities_raster_int(:,ind_n),ground_truth);
        
        %plot
        if pintar_all && ind_n<16
            %plot traces
            subplot(4,4,ind_n)
            hold on
            if noise_photon_level==0 && neuropil_level==0
                plot(frame_times_scan_path(:,ind_n)/1000,activities_scan_path(:,ind_n),'+-c','linewidth',1)
                plot(timing_raster(:,ind_n)/1000,activities_raster(:,ind_n),'+-r','linewidth',2)
                plot(1/1000:1/1000:data_toy.duration_ms/1000,ground_truth,'k','linewidth',2)
            else
                plot(frame_times_scan_path(:,ind_n)/1000,activities_scan_path(:,ind_n)-mean(activities_scan_path(:,ind_n)),'+-c','linewidth',1)
                plot(timing_raster(:,ind_n)/1000,activities_raster(:,ind_n)-mean(activities_raster(:,ind_n)),'+-r','linewidth',2)
                plot(1/1000:1/1000:data_toy.duration_ms/1000,ground_truth-mean(ground_truth),'k','linewidth',2)
            end
            title(['scan corr: ' num2str(scan_corr(ind_n),2) ' / raster corr: ' num2str(raster_corr(ind_n),2)])
        end
        %plot example trace for summary figure
        if pintar && ind_n==1 && mat_surrounding(ind_surr)==path_panel_num_pixels
            hold(axes_stats(1),'on')
            plot(axes_stats(1),frame_times_scan_path(:,ind_n)/1000,activities_scan_path(:,ind_n)-mean(activities_scan_path(:,ind_n)),'+-c','linewidth',1)
            plot(axes_stats(1),timing_raster(:,ind_n)/1000,activities_raster(:,ind_n)-mean(activities_raster(:,ind_n)),'+-r','linewidth',2)
            plot(axes_stats(1),1/1000:1/1000:data_toy.duration_ms/1000,ground_truth-mean(ground_truth),'k','linewidth',2)
            plot(axes_stats(1),1/1000:1/1000:data_toy.duration_ms/1000,...
                ground_truth_noise_and_neuropil-mean(ground_truth_noise_and_neuropil),'color',[.5 .5 .5],'linewidth',2)
            xlabel(axes_stats(1),'time (s)')
            ylabel(axes_stats(1),'fluorescence (a.u.)')
            title(axes_stats(1),['scan corr: ' num2str(scan_corr(ind_n),2) ' (' num2str(frame_rates(ind_surr),2)...
                'Hz) / raster corr: ' num2str(raster_corr(ind_n),2)])
            legend(axes_stats(1),'smart-line scan','raster scan','ground truth','ground truth with neuropil')
            set(axes_stats(1),'fontSize',20)
        end
    end
    if pintar_all
        %plot neuropil activity
        subplot(4,4,16)
        hold on
        ground_truth = mat_fluorescence(end,:);
        ground_truth = ground_truth' - mean(ground_truth);
        plot(1/1000:1/1000:data_toy.duration_ms/1000,ground_truth,'k','linewidth',2)
        title('neuropil')
    end
    %save correlation
    mean_performance_raster(:,ind_surr) = raster_corr;
    mean_performance_scan(:,ind_surr) = scan_corr;
    
    
    
    %measure SNRChen (as in Chen et al. 2013)
    for ind_n=1:num_neurons
        %timings
        raster_time = timing_raster(:,ind_n);
        scan_time = frame_times_scan_path(:,ind_n);
        pintar_SNRChen = ind_n==1 && mod(ind_surr,2)==0 && pintar_all;
        %in case there is no stim (this part of the code will probably be removed at some point)
        if isempty(stim_times)
            if pintar_SNRChen
                set(0,'currentFigure',h13)
            end
            %get the trace
            trace = activities_raster(:,ind_n);
            stim_time = 0;
            %get SNRCHen values for raster
            [SNRChen_raster(ind_n,ind_surr),SNRChen_raster_max(ind_n,ind_surr),SNRChen_raster_std(ind_n,ind_surr)] =...
                SNRChen_calculation(trace,25,10000, stim_time, pintar_SNRChen*panel_counter,raster_time);
            
            if pintar_SNRChen
                set(0,'currentFigure',h14)
            end
            %the SNRChen values shown in the panels of this figure corresponds
            %only to the particular trial displayed and will normally not macth
            %the final SNRChen of the neuron
            %get the trace
            trace = activities_scan_path(:,ind_n);
            stim_time = 0;
            [SNRChen_scan(ind_n,ind_surr),SNRChen_max(ind_n,ind_surr),SNRChen_std(ind_n,ind_surr)] = SNRChen_calculation(trace,25,10000,stim_time, pintar_SNRChen*panel_counter,scan_time);
        else
            %go over all stimuli
            for ind_stim=1:numel(stim_times)
                %raster
                if pintar_SNRChen
                    set(0,'currentFigure',h13)
                end
                
                %get trace and time points corresponding to the current stimulus
                trace = activities_raster(:,ind_n);
                index = raster_time>stim_times(ind_stim)-baseline_used_for_exps & raster_time<=stim_times(ind_stim)+response_used_for_exps;
                trace = trace(index);
                stim_time = find(raster_time>stim_times(ind_stim),1,'first')-find(raster_time>stim_times(ind_stim)-baseline_used_for_exps,1,'first');
                [SNRChen_raster_aux,SNRChen_raster_max_aux,SNRChen_raster_std_aux] = SNRChen_calculation(trace,25,10000, stim_time, pintar_SNRChen*(ind_stim==1)*panel_counter,raster_time(index));
                SNRChen_raster(ind_n,ind_surr) = SNRChen_raster(ind_n,ind_surr) + SNRChen_raster_aux;
                SNRChen_raster_max(ind_n,ind_surr) = SNRChen_raster_max(ind_n,ind_surr) + SNRChen_raster_max_aux;
                SNRChen_raster_std(ind_n,ind_surr) = SNRChen_raster_std(ind_n,ind_surr) + SNRChen_raster_std_aux;
                %lineScan
                if pintar_SNRChen
                    set(0,'currentFigure',h14)
                end
                %get trace and time points corresponding to the current stimulus
                trace = activities_scan_path(:,ind_n);
                index = scan_time>stim_times(ind_stim)-baseline_used_for_exps & scan_time<=stim_times(ind_stim)+response_used_for_exps;
                trace = trace(index);
                stim_time = find(scan_time>stim_times(ind_stim),1,'first')-find(scan_time>stim_times(ind_stim)-baseline_used_for_exps,1,'first');
                [SNRChen_scan_aux,SNRChen_max_aux,SNRChen_std_aux] = SNRChen_calculation(trace,25,10000,stim_time, pintar_SNRChen*(ind_stim==1)*panel_counter,scan_time(index));
                SNRChen_scan(ind_n,ind_surr) = SNRChen_scan(ind_n,ind_surr) + SNRChen_scan_aux;
                SNRChen_max(ind_n,ind_surr) = SNRChen_max(ind_n,ind_surr) + SNRChen_max_aux;
                SNRChen_std(ind_n,ind_surr) = SNRChen_std(ind_n,ind_surr) + SNRChen_std_aux;
                %ground truth SNRChen
                %get trace and time points corresponding to the current stimulus
                trace = activities_ground_truth_noise_and_neuropil(:,ind_n);
                timing = 1:numel(ground_truth);
                index = timing>stim_times(ind_stim)-baseline_used_for_exps & timing<=stim_times(ind_stim)+response_used_for_exps;
                trace = trace(index);
                stim_time = find(timing>stim_times(ind_stim),1,'first')-find(timing>stim_times(ind_stim)-baseline_used_for_exps,1,'first');
                [SNRChen_scan_aux,~,~] = SNRChen_calculation(trace,25,10000,stim_time, 0,timing(index));
                SNRChen_gt(ind_n,ind_surr) = SNRChen_gt(ind_n,ind_surr) + SNRChen_scan_aux;
            end
            %average
            SNRChen_raster(ind_n,ind_surr) = SNRChen_raster(ind_n,ind_surr)/numel(stim_times);
            SNRChen_raster_max(ind_n,ind_surr) = SNRChen_raster_max(ind_n,ind_surr)/numel(stim_times);
            SNRChen_raster_std(ind_n,ind_surr) = SNRChen_raster_std(ind_n,ind_surr)/numel(stim_times);
            SNRChen_scan(ind_n,ind_surr) = SNRChen_scan(ind_n,ind_surr)/numel(stim_times);
            SNRChen_max(ind_n,ind_surr) = SNRChen_max(ind_n,ind_surr)/numel(stim_times);
            SNRChen_std(ind_n,ind_surr) = SNRChen_std(ind_n,ind_surr)/numel(stim_times);
            SNRChen_gt(ind_n,ind_surr) = SNRChen_gt(ind_n,ind_surr)/numel(stim_times);
        end
        panel_counter = panel_counter + pintar_SNRChen;
    end
    
    
    %measure SNRMarco (from single spikes)
    for ind_n=1:num_neurons
        %get selected spikes
        spiketimes = ground_truth_spiketimes{ind_n};
        aux = diff([0 spiketimes duration_ms]);
        select_spikes = aux(1:end-1)>baseline_window & aux(2:end)>peak_window;
        spiketimes(~select_spikes) = [];
        %get trace and scan times
        trace = activities_scan_path(:,ind_n);
        scan_time = frame_times_scan_path(:,ind_n);
        if ~isempty(spiketimes)
            [SNRMarco_scan(ind_n,ind_surr),SNRMarco_max(ind_n,ind_surr),SNRMarco_std(ind_n,ind_surr)] =...
                SNRMarco_calculation(trace, spiketimes, scan_time, baseline_window, peak_window);
        end
    end
    
    if pintar && ind_surr==numel(mat_surrounding)
        hold(axes_stats(2),'on')
        %plot results regarding correlation VS frame rate
        mean_raster = mean(nanmax(mean_performance_raster,[],2));
        plt_1 = plot(axes_stats(2),[frame_rates(1),frame_rates(ind_surr)],mean_raster*ones(1,2),'--r','lineWidth',2);
        plt_2 = plot(axes_stats(2),frame_rates(1:ind_surr),mean_performance_scan(:,1:ind_surr)','color',[.5 .5 1]);
        plt_3 = plot(axes_stats(2),frame_rates(1:ind_surr),nanmean(mean_performance_scan(:,1:ind_surr),1),'b','lineWidth',2);
        xlabel(axes_stats(2),'sampling rate (Hz)')
        ylabel(axes_stats(2),'correlation with ground truth')
        legend(axes_stats(2),[plt_1,plt_2(1),plt_3],{'mean corr. raster scan','mean corr. smart-line scan','ind. n. smart-line scan'},'Location','southwest')
        set(axes_stats(2),'fontSize',20)
        axis tight
    end
    
    if pintar_all
        if ind_surr==numel(mat_surrounding)
            %plot SNRChen (maybe this part of the code could be removed)
            set(0,'currentFigure',h12)
            subplot(1,3,1)
            hold on
            plot(frame_rates(1:ind_surr),SNRChen_scan(2:end,1:ind_surr)','color',[.5 .5 1])
            plot(frame_rates(1:ind_surr),nanmean(SNRChen_scan(:,1:ind_surr),1),'b','lineWidth',2)
            plot(frame_rates(1:ind_surr),nanmean(SNRChen_raster(:,1:ind_surr),1),'r','lineWidth',2)
            plot(frame_rates(1:ind_surr),SNRChen_scan(1,1:ind_surr)','color',[1 0 0])
            xlabel(gca,'sampling rate (Hz)')
            ylabel(gca,'SNRChenChen')
            set(gca,'fontSize',20)
            title('SNRChenChen')
            axis square
            hold off
            
            
            subplot(1,3,2)
            hold on
            plot(frame_rates(1:ind_surr),SNRChen_max(2:end,1:ind_surr)','color',[.5 .5 1])
            plot(frame_rates(1:ind_surr),nanmean(SNRChen_max(:,1:ind_surr),1),'b','lineWidth',2)
            plot(frame_rates(1:ind_surr),SNRChen_max(1,1:ind_surr)','color',[1 0 0])
            xlabel(gca,'sampling rate (Hz)')
            ylabel(gca,'dynamic range')
            set(gca,'fontSize',20)
            title('dynamic range')
            axis square
            hold off
            
            subplot(1,3,3)
            hold on
            plot(frame_rates(1:ind_surr),SNRChen_std(2:end,1:ind_surr)','color',[.5 .5 1])
            plot(frame_rates(1:ind_surr),nanmean(SNRChen_std(:,1:ind_surr),1),'b','lineWidth',2)
            plot(frame_rates(1:ind_surr),SNRChen_std(1,1:ind_surr)','color',[1 0 0])
            xlabel(gca,'sampling rate (Hz)')
            ylabel(gca,'fluorescence std')
            set(gca,'fontSize',20)
            title('fluorescence std')
            axis square
            hold off
            
            
            %plot individual neurons performances
            set(0,'currentFigure',h7)
            for ind_n=1:num_neurons
                subplot(ceil(sqrt(num_neurons)),ceil(sqrt(num_neurons)),ind_n)
                if ind_surr==numel(mat_surrounding)
                    plot([frame_rates(1),frame_rates(ind_surr)],max(mean_performance_raster(ind_n,:))*ones(1,2),'--','color','r')
                end
                hold on
                plot(frame_rates(1:ind_surr),mean_performance_scan(ind_n,1:ind_surr)','color','b')
                hold off
            end
        end
        close(h3)
        close(h4)
        close(h6)
    end
    save([folder_name 'data' num2str(mat_surrounding(ind_surr))],'whole_path','activities_scan_path',...
        'frame_times_scan_path','tags','activities_raster','timing_raster','mean_performance_scan',...
        'mean_performance_raster','activities_ground_truth','activities_ground_truth_noise_and_neuropil',...
        'effective_pixels_per_neuron','effective_num_pixels_per_neuron')
end
%save figures
if pintar
    for ind=1:4
        set(axes_rasters(ind),'xtick',[])
        set(axes_rasters(ind),'ytick',[])
    end
    saveas(h_paper,[folder_name 'summary figure'],'fig')
    print(h_paper,[folder_name 'summary figure'],'-dpng')
    saveas(h_paper,[folder_name 'summary figure'],'epsc')
end
if pintar_all
    
    %%%%%%
    saveas(h12,[folder_name '_SNRChen'],'fig')
    print(h12,[folder_name '_SNRChen'],'-dpng')
    saveas(h12,[folder_name '_SNRChen'],'epsc')
    %%%%%%
    saveas(h13,[folder_name '_SNRChen_example_neuron_raster'],'fig')
    print(h13,[folder_name '_SNRChen_example_neuron_raster'],'-dpng')
    saveas(h13,[folder_name '_SNRChen_example_neuron_raster'],'epsc')
    %%%%%%
    saveas(h14,[folder_name '_SNRChen_example_neuron'],'fig')
    print(h14,[folder_name '_SNRChen_example_neuron'],'-dpng')
    saveas(h14,[folder_name '_SNRChen_example_neuron'],'epsc')
    
    saveas(h8,[folder_name 'movie'],'fig')
    saveas(h8,[folder_name 'movie'],'png')
    saveas(h7,[folder_name 'performances neurons different N'],'fig')
    print(h7,[folder_name 'performances neurons different N'],'-dpng')
    
end

%put results on the data_toy struct
data_toy.frame_rates = frame_rates;
data_toy.SNRChen_scan = SNRChen_scan;
data_toy.SNRChen_max = SNRChen_max;
data_toy.SNRChen_std = SNRChen_std;
data_toy.SNRChen_raster_scan = SNRChen_raster;
data_toy.SNRChen_raster_max = SNRChen_raster_max;
data_toy.SNRChen_raster_std = SNRChen_raster_std;
data_toy.SNRChen_gt = SNRChen_gt;
data_toy.SNRMarco_scan = SNRMarco_scan;
data_toy.SNRMarco_max = SNRMarco_max;
data_toy.SNRMarco_std = SNRMarco_std;

data_toy.performance_raster = mean_performance_raster;
data_toy.performance_scan = mean_performance_scan;
save([folder_name name],'data_toy')
disp('sampling rates')
disp(data_toy.frame_rates)
disp('mean transient amplitude')
disp(nanmean(SNRMarco_max,1))

close all
end

function [trace,fr_response,fr_baseline] = circuit_II(t_max,spike_rate, stim_times, baseline_rate)
%this will create the firing rate trace and return a spike train based on it and a poisson process
if isempty(stim_times)
    spike_rate_mat = spike_rate*ones(1,t_max)/1000;
else
    response_duration = stim_times(2)-2*stim_times(1);
    response = 0:response_duration-1;
    baseline_spike_rate = baseline_rate/1000;
    peak_response_spike_rate = spike_rate/1000;
    response = (peak_response_spike_rate-baseline_spike_rate)*exp(-response/10) + baseline_spike_rate;%response*(baseline_spike_rate-peak_response_spike_rate)/response(end) + peak_response_spike_rate;
    spike_rate_mat = zeros(1,stim_times(end)+response_duration);
    trial_duration = (stim_times(1)+response_duration);
    for ind=1:numel(stim_times)
        spike_rate_mat((ind-1)*trial_duration+1:ind*trial_duration) = [baseline_spike_rate*ones(1,stim_times(1)) response];
    end
    spike_rate_mat = [spike_rate_mat ones(1,t_max-numel(spike_rate_mat))*baseline_spike_rate];
end

trace = poissrnd(spike_rate_mat);
trace(trace>0) = 1;
fr_response = 1000*trace(stim_times(1)+1:stim_times(1)+50);
fr_baseline = 1000*nnz(trace(1:stim_times(1)))/stim_times(1);
%save figure with firing rate mat
%{
h = figure('units','normalized','position',[0.1 0.1 0.8 0.4]);
hold on
baseline_duration = 50;
plot(-baseline_duration:100,1000*spike_rate_mat(stim_times(1)-baseline_duration:stim_times(1)+100),'lineWidth',2)
%plot((baseline_duration-20)*ones(1,2),[0 25],'--','lineWidth',2,'color',[.6 .6 .6])
plot(zeros(1,2),[0 25],'--','lineWidth',2,'color',[.6 .6 .6])
plot(50*ones(1,2),[0 25],'--','lineWidth',2,'color',[.6 .6 .6])
title(['Mean firing rate (50ms post-stim): ' num2str((1000*mean(spike_rate_mat(stim_times(1)+1:stim_times(1)+50))),2) 'Hz'])
xlabel('time (ms)')
ylabel('firing rate (Hz)')
set(gca,'fontSize',20)
xlim([-50 100])
%}
end

function control = control_vecindad(ind1,ind2,radio,Mos)
%this is just to control that a neuron is not too close to any of the rest
dist = sqrt((ind1-Mos(1,:)).^2+(ind2-Mos(2,:)).^2);
if min(dist)<=radio
    control = 1;
else
    control = 0;
end

end

function mos_random = mosaico_random_induccion(tam,sigma_tam_neuronas,num_neurons)
%returns a matrix with the the position of each neuron
vecindad = 8*sigma_tam_neuronas;
num = num_neurons;
stop = false;
mos_random = zeros(2,num);
for s = 1:num
    i = unidrnd(tam);
    j = unidrnd(tam);
    counter = 0;
    while ((control_vecindad(i,j,vecindad,mos_random)) || i<vecindad || i>tam-vecindad || j<vecindad || j>tam-vecindad) && ~stop
        i = unidrnd(tam);
        j = unidrnd(tam);
        counter = counter + 1;
        stop = counter>1000;
    end
    if ~stop
        mos_random(:,s)= [i,j]';
    else
        break
    end
end
if stop
    final_num_neurons = 2*floor((s-1)/2);
else
    final_num_neurons = 2*floor(s/2);
end
mos_random = mos_random(:,1:final_num_neurons);
assert(nnz(mos_random==0)==0)
disp(['population of ' num2str(final_num_neurons) ' neurons (expected: ' num2str(num_neurons) ')'])
end

function [r,pixels_patch] = neuropil_activity(duration_ms)
close all
folder = '/home/manuel/Desktop/juxta and imaging/results_summary/';
load([folder 'final_results.mat'])
% figure
% imagesc(all_cov_mat)
r = mvnrnd(nanmean(all_mean_fluorescences)*zeros(1,size(all_cov_mat,1)),all_cov_mat,duration_ms);%note that this has zero mean!!
pixels_patch = pixels_patch - round(neuropil_patch_size/2); %#ok<NODEF>
end

function kernel = smoothing_kernel(index,smoothing_factor)
%creates a kernel that will be used to smooth the neurons/neuropil shape
mu = [0 0];
Sigma = [smoothing_factor 0; 0 smoothing_factor];
[X1,X2] = meshgrid(index,index);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
kernel = reshape(F,length(index),length(index));
kernel = kernel/sum(kernel(:));
end

function [movie_sample,movie_sample_time] = movie_sampling_raster(movie_final,frame_times,dwellTime,tam,noise_photon_level)
frame_times = 1000*frame_times;%frame times in ms
dwellTime = 1000*dwellTime;%dwell time in ms
contador = 1;
movie_sample = zeros(numel(frame_times),tam^2);
movie_sample_time = zeros(numel(frame_times),size(movie_final,2));
average_movie = mean(movie_final,1);
for ind_px_1=1:tam
    for ind_px_2=1:tam
        indx = sub2ind([tam tam],ind_px_1,ind_px_2);
        %get activity
        aux = movie_final(floor(frame_times+(ind_px_1-1)*tam*dwellTime+(ind_px_2-1)*dwellTime)+1,indx);
        
        %add noise
        random_noise = randn(size(aux));
        %sampling noise depends on fluorescence level (see Pnevmatikakis et
        %al.). The +1 is just to avoid negative values
        std_noise = noise_photon_level*average_movie(indx);
        random_noise = random_noise.*std_noise+1;
        
        %save timing
        movie_sample_time(:,indx) = floor(frame_times+(ind_px_1-1)*tam*dwellTime+(ind_px_2-1)*dwellTime)+1;
        movie_sample(:,contador) = aux + random_noise(:);
        contador = contador + 1;
    end
end
%remove negative values
movie_sample(movie_sample<0) = 0;
end

function imagen = find_strongest_event(video,tam,all_pixels_on_square,average_window_choose_roi_current,duration)
%find strongest event in an area determined by all_pixels_on_square
global_activity = sum(video(:,all_pixels_on_square),2);
[~,event] = max(global_activity);
current_stack = max([event-average_window_choose_roi_current,1]):...
    min([event+average_window_choose_roi_current,duration]);
imagen = mean(video(current_stack,:),1);

imagen = reshape(imagen,tam,[]);
end

function [scan_path, tags, num_selected_pixels, signal_to_noise_ratio_mat] =...
    get_SNRGUI_for_different_num_of_pixels(imagen,all_pixels_on_square_2d, all_pixels_on_square, video, ind_n, scan_path, num_pixels_th, counter, tags, pintar)
%get pixels associated with a given neuron by measuing the SNR for increasing number of pixels (ordered by brightness)
num_figures = 1;
num_steps = 150;
num_steps = min(num_steps,numel(all_pixels_on_square));
fluorescences = imagen(all_pixels_on_square);
[~,index_fl] = sort(fluorescences,'descend');
all_pixels_on_square_sorted = all_pixels_on_square(index_fl);
all_pixels_on_square_2D_sorted = all_pixels_on_square_2d(index_fl,:);
signal_to_noise_ratio_mat = nan(1,150);
if pintar && ind_n<num_figures
    figure('units','normalized','position',[.05 .05 .9 .9]);
end
mat_num_pixels = 5:num_steps;%I do not allow selecting less than 5 pixels 
for ind_px=mat_num_pixels
    %get trace
    trace = mean(video(:,all_pixels_on_square_sorted(1:ind_px)),2);
    %compute SNR (from max and std)
    [maximo,argmaximo] = max(trace);
    lower_trace_th_aux = prctile(trace,25);
    lower_trace_activity = trace(trace<=lower_trace_th_aux);
    if std(lower_trace_activity)==0
        lower_standard_deviation = inf;
    else
        lower_standard_deviation = std(lower_trace_activity);
    end
    signal_to_noise_ratio_mat(ind_px) = (max(trace)-mean(lower_trace_activity))/lower_standard_deviation;
    if pintar && mod(ind_px,2)==0 && ind_n<num_figures
        subplot(10,10,ind_px/2)
        hold on
        t = 1:numel(trace);
        plot(t,trace)
        plot(t(trace<=lower_trace_th_aux),lower_trace_activity,'*c')
        plot(argmaximo,maximo,'+m')
        set(gca,'xtick',[],'ytick',[])
        title(['SNR=' num2str(signal_to_noise_ratio_mat(ind_px))])
    end
end
%define the threshold used to select the number of pixels
signal_to_noise_ratio_mat_norm = (signal_to_noise_ratio_mat-nanmin(signal_to_noise_ratio_mat))/(nanmax(signal_to_noise_ratio_mat)-nanmin(signal_to_noise_ratio_mat));
SNR_threshold = num_pixels_th;%nanmax(signal_to_noise_ratio_mat)-(nanmax(signal_to_noise_ratio_mat)-nanmin(signal_to_noise_ratio_mat))/10;
argmaximo = find(signal_to_noise_ratio_mat_norm>SNR_threshold,1,'first');
maximo = signal_to_noise_ratio_mat_norm(argmaximo);

%get selected pixels and put them on the scan path
selected_pixels = all_pixels_on_square_2D_sorted(1:argmaximo,:);
num_selected_pixels = size(selected_pixels,1);
scan_path(:,counter+1:counter+num_selected_pixels) =...
    selected_pixels';
tags(counter+1:counter+num_selected_pixels) = ind_n;
assert(num_selected_pixels>0)
if pintar && ind_n<num_figures
    figure('units','normalized','position',[.1 .1 .8 .8]);
    subplot(2,6,1:6)
    hold on
    plot(1:numel(signal_to_noise_ratio_mat_norm),signal_to_noise_ratio_mat_norm,'lineWidth',2)
    plot([1,numel(signal_to_noise_ratio_mat_norm)],SNR_threshold*ones(1,2),'--c')
    legend('original','threshold')
    for ind=1:5
        subplot(2,6,6+ind)
        imagesc(imagen)
        hold on
        plot(all_pixels_on_square_2D_sorted(1:ind*floor(num_steps/6),2),all_pixels_on_square_2D_sorted(1:ind*floor(num_steps/6),1),'+r')
        title(['num pixels=' num2str(ind*round(num_steps/6))])
        xlim([min(all_pixels_on_square_2D_sorted(:,2))-1,max(all_pixels_on_square_2D_sorted(:,2))+1])
        ylim([min(all_pixels_on_square_2D_sorted(:,1))-1,max(all_pixels_on_square_2D_sorted(:,1))+1])
        subplot(2,6,1:6)
        plot(ind*round(num_steps/6)*ones(1,2),[nanmin(signal_to_noise_ratio_mat_norm) nanmax(signal_to_noise_ratio_mat_norm)],'--k')
    end
    plot(argmaximo,maximo,'*r','markerSize',12)
    
    subplot(2,6,12)
    imagesc(imagen)
    hold on
    set(gca,'xcolor','r','ycolor','r','linewidth',2)
    plot(all_pixels_on_square_2D_sorted(1:argmaximo,2),all_pixels_on_square_2D_sorted(1:argmaximo,1),'+r')
    title(['num pixels=' num2str(num_selected_pixels)])
    xlim([min(all_pixels_on_square_2D_sorted(:,2))-1,max(all_pixels_on_square_2D_sorted(:,2))+1])
    ylim([min(all_pixels_on_square_2D_sorted(:,1))-1,max(all_pixels_on_square_2D_sorted(:,1))+1])
end
end

function [scan_path, tags_neurons, num_selected_pixels] =...
    get_SNRGUI_for_different_scanLines(transient, all_pixels_on_square_2d, ind_n, scan_path, counter, tags_neurons, tam)
%get pixels associated with a given neuron by measuing the SNR for increasing number of pixels (ordered by brightness)
x0_mat = sort(unique(all_pixels_on_square_2d(:,1)));
y0_mat = unique(all_pixels_on_square_2d(:,2));
min_x = min(x0_mat);
max_x = max(x0_mat);
min_y = min(y0_mat);
max_y = max(y0_mat);
x0_mat = x0_mat(2:end-1);%remove the e
%add all pixels in the borders
aux = [x0_mat,min_y*ones(numel(x0_mat),1)];
aux = unique(aux,'rows');
all_pixels = aux;
tags = zeros(1,size(aux,1));
aux = [x0_mat,max_y*ones(numel(x0_mat),1)];
aux = unique(aux,'rows');
all_pixels = [all_pixels; aux];
tags = [tags zeros(1,size(aux,1))+1];

aux = [min_x*ones(numel(y0_mat),1),y0_mat];
aux = unique(aux,'rows');
all_pixels = [all_pixels; aux];
tags = [tags zeros(1,size(aux,1))+2];

aux = [max_x*ones(numel(y0_mat),1),y0_mat];
aux = unique(aux,'rows');
all_pixels = [all_pixels; aux];
tags = [tags zeros(1,size(aux,1))+3];
final_fluorescence_captured = 0;
for ind_border=0:3
    pixels1 = all_pixels(tags==ind_border,:);
    for ind_borderII=ind_border+1:3
        pixels2 = all_pixels(tags==ind_borderII,:);
        
        for ind_p1=1:size(pixels1)
            point1 = pixels1(ind_p1,:);
            for ind_p2=1:size(pixels2)
                point2 = pixels2(ind_p2,:);
                if abs(point2(1)-point1(1))>abs(point2(2)-point1(2))
                    x_mat = min([point1(1),point2(1)]):0.1:max([point1(1),point2(1)]);
                    y_mat = (point2(2)-point1(2))*(x_mat-point1(1))/(point2(1)-point1(1))+point1(2);
                    line = round([x_mat;y_mat])';
                else
                    x_mat = min([point1(2),point2(2)]):0.1:max([point1(2),point2(2)]);
                    y_mat = (point2(1)-point1(1))*(x_mat-point1(2))/(point2(2)-point1(2))+point1(1);
                    line = round([y_mat;x_mat])';
                end
                
                line = unique(line,'rows','stable');
                line_1D = sub2ind([tam tam],line(:,1),line(:,2));
                fluorescence = transient(line_1D);
                line(fluorescence<prctile(fluorescence,5),:) = [];
                fluorescence(fluorescence<prctile(fluorescence,5)) = [];
                fluorescence_captured = sum(fluorescence)/numel(fluorescence);
                
                %get the pixels in the scan path intersecting with pixels_in
                %assert(~isempty(intersect(line,all_pixels_on_square_2d,'rows')))
                
                if final_fluorescence_captured<fluorescence_captured && ~isempty(intersect(line,all_pixels_on_square_2d,'rows'))
                    final_fluorescence_captured = fluorescence_captured;
                    selected_pixels = line;
                end
                %{
                 figure
                imagesc(transient)
                hold on
                plot(selected_pixels(:,1),selected_pixels(:,2),'r')
                plot(line(:,1),line(:,2))
                ylim([min_x-5 max_x+5])
                xlim([min_y-5 max_y+5])
                title([num2str(fluorescence_captured) ' (best = ' num2str(signal_to_noise_ratio) ')'])
                close all
                %}
            end
        end
    end
end

num_selected_pixels = size(selected_pixels,1);
%get selected pixels and put them on the scan path
scan_path(:,counter+1:counter+num_selected_pixels) =...
    selected_pixels';
tags_neurons(counter+1:counter+num_selected_pixels) = ind_n;
assert(num_selected_pixels>0)

%{
h = figure;
imagesc(transient)
hold on
plot(selected_pixels(:,1),selected_pixels(:,2),'r','lineWidth',2)
ylim([min_x-5 max_x+5])
xlim([min_y-5 max_y+5])
close(h)
%}
end

function [activity,frame_times] = get_raster_activities(movie,scan_path,tags,ind,tam,raster_sampling_time)
%returns activities for raster from the pixels specified in scan_path/tags.
%If sampling_time is provided it returns also the frame times.
indx_tags = find(tags==ind);
indx = sub2ind([tam tam],scan_path(1,indx_tags),scan_path(2,indx_tags));
activity = mean(movie(:,indx),2);

if nargin>5
    %the dimensions are swapped for the raster_sampling_time matrix
    indx = sub2ind([tam tam],scan_path(2,indx_tags),scan_path(1,indx_tags));
    frame_times = round(mean(raster_sampling_time(:,indx),2))';
    assert(nnz(isnan(frame_times))==0)
end
end

function whole_path = include_all_pixels(scan_path,tam,pintar)
%add in between pixels to a path (the prairie laser cannot jump)
num_pixels = size(scan_path,2);
whole_path = zeros(2,tam^2);
whole_path(:,1) = scan_path(:,1);
last_point = whole_path(:,1);
contador = 1;
for ind_px=2:num_pixels
    x = last_point;
    y = scan_path(:,ind_px);
    if pintar
        plot(y(2),y(1),'+y')
    end
    if abs(x(1)-y(1))>abs(x(2)-y(2))
        t = x(1):sign(y(1)-x(1))*0.1:y(1);
        line_aux = (y(2)-x(2))*(t-x(1))/(y(1)-x(1))+x(2);
    else
        line_aux = x(2):sign(y(2)-x(2))*0.1:y(2);
        t = (y(1)-x(1))*(line_aux-x(2))/(y(2)-x(2))+x(1);
    end
    inter_pixel_path = unique(round([t;line_aux])','rows','stable')';
    whole_path(:,contador+1:contador+size(inter_pixel_path,2)+1) = [inter_pixel_path y];
    contador = contador + size(inter_pixel_path,2) + 1;
    last_point = y;
end
whole_path = whole_path(:,1:contador);
%find all moments in which the laser would stop in one pixel
diff_ = sum(abs(diff(whole_path,[],2)),1);
whole_path = whole_path(:,diff_~=0);

end

function path_all = group_pixels_and_calculate_shortest_path_lineScan(path,tags)
%this function groups calculate each pixel centroid and call shortest_path
%to get the shortest path it can find. Then it builds the whole path by
%applying a greedy search to each neuron's pixels.
% close all
% tic
tags_unique = unique(tags);
nodes = zeros(numel(tags_unique),2);
for ind_n=1:numel(tags_unique)
    nodes(ind_n,:) = mean(path(tags==tags_unique(ind_n),:),1);
end


userConfig = struct('xy',nodes,'showProg',false,'showResult',false,'showWaitbar',false);
resultStruct = tsp_ga(userConfig);
path_neurons = resultStruct.xy(resultStruct.optRoute,:);% [path_neurons,distancia] = shortest_path(nodes);

%once we get the path going over all neurons, we put all pixels from all neurons together
tag_aux = tags_unique(sum(abs(nodes-repmat(path_neurons(1,:),numel(tags_unique),1)),2)==0);
nodes_group = path(tags==tag_aux,:);

distance1 = sum(sqrt(sum((nodes_group(1,:)-path).^2,2)));
distance2 = sum(sqrt(sum((nodes_group(end,:)-path).^2,2)));
if distance1>distance2
    path_all = nodes_group;
    last_point = nodes_group(end,:);
else
    path_all = nodes_group(end:-1:1,:);
    last_point = nodes_group(1,:);
end

for ind_tag=1:numel(tags_unique)
    
    tag_aux = tags_unique(sum(abs(nodes-repmat(path_neurons(ind_tag,:),numel(tags_unique),1)),2)==0);
    nodes_group = path(tags==tag_aux,:);
    
    distance1 = sum(sqrt(sum((nodes_group(1,:)-last_point).^2,2)));
    distance2 = sum(sqrt(sum((nodes_group(end,:)-last_point).^2,2)));
    if distance1<distance2
        path_all = [path_all;nodes_group]; %#ok<AGROW>;
        last_point = nodes_group(end,:);
    else
        path_all = [path_all;nodes_group(end:-1:1,:)]; %#ok<AGROW>;
        last_point = nodes_group(1,:);
    end
end
end

function [movie_sample,movie_sample_only_neurons,movie_sample_time] =...
    scan_path_movie(whole_path,dwellTime,data_toy,noise_photon_level,movie_ruido,movie)
%scan movie with the SmartLine scan
num_pixels_total = size(whole_path,2);
framePeriod = dwellTime*num_pixels_total;
frame_times =  0:framePeriod:data_toy.duration_ms/1000-framePeriod;
frame_times = 1000*frame_times;%frame times in ms
dwellTime = 1000*dwellTime;
movie_sample = zeros(numel(frame_times),num_pixels_total);
movie_sample_only_neurons = zeros(numel(frame_times),num_pixels_total);
movie_sample_time = zeros(numel(frame_times),num_pixels_total);
average_movie = mean(movie_ruido,1);
for ind_px=1:num_pixels_total
    indx = sub2ind([data_toy.tam data_toy.tam],whole_path(2,ind_px),whole_path(1,ind_px));
    aux = movie_ruido(floor(frame_times+(ind_px-1)*dwellTime)+1,indx);
    
    %add noise
    random_noise = randn(size(aux));
    %sampling noise depends on fluorescence level (see Pnevmatikakis et
    %al.). The +1 is just to avoid negative values
    std_noise = noise_photon_level*average_movie(indx);
    random_noise = random_noise.*std_noise+1;
    
    movie_sample(:,ind_px) = aux + random_noise(:);
    
    aux = movie(floor(frame_times+(ind_px-1)*dwellTime)+1,indx);
    movie_sample_only_neurons(:,ind_px) = aux;
    
    %save alos the times
    movie_sample_time(:,ind_px) = floor(frame_times+(ind_px-1)*dwellTime)+1;
end
movie_sample(movie_sample<0) = 0;
end

function [activities, frame_times_scan,tags, spike_times, effective_pixels_per_neuron, effective_num_pixels_per_neuron] =...
    get_rois_and_activities(movie_sample, pixels_per_neuron, num_neurons, scan_path, movie_sample_time)
%segment the scan movie and get the activities
num_pixels_total = size(scan_path,2);

activities = zeros(size(movie_sample,1),num_neurons);
frame_times_scan = zeros(size(movie_sample,1),num_neurons);
spike_times = cell(1,num_neurons);
tags = zeros(1,num_pixels_total);
effective_pixels_per_neuron = cell(1,num_neurons);
effective_num_pixels_per_neuron = zeros(1,num_neurons);
for ind_n=1:num_neurons
    %get the pixels close to the current neuron
    pixels_in = pixels_per_neuron{ind_n};%select_pixels(positions(:,ind_n), tam, reference, sigma_tam_neuronas, sigma_factor);
    %get the pixels in the scan path intersecting with pixels_in
    scan_path_aux = scan_path';
    %intersect does not return ALL pixels in the path that pass through the
    %ones in pixels_in, and so here I am making sure that we take them all
    [~,indx_aux] = intersect(scan_path_aux,pixels_in,'rows');
    indx = indx_aux;
    while ~isempty(indx_aux)
        scan_path_aux(indx_aux,:) = nan;
        [~,indx_aux] = intersect(scan_path_aux,pixels_in,'rows');
        indx = [indx; indx_aux]; %#ok<AGROW>
    end
    
    assert(~isempty(indx))
    effective_pixels_per_neuron{ind_n} = scan_path(:,indx);
    effective_num_pixels_per_neuron(ind_n) = numel(indx);
    activities(:,ind_n) = mean(movie_sample(:,indx),2);
    
    frame_times_scan(:,ind_n) = round(mean(movie_sample_time(:,indx),2));
    
    tags(indx) = ind_n;
    
    
end
end

function  activities_int = interpolate_activity(activities,frame_times,duration)
%interpolate fluorescence activity so it can be compared to ground truth
activities_int = zeros(duration,size(activities,2));

for ind_n=1:size(activities,2)
    activity = activities(:,ind_n);
    timing = frame_times(:,ind_n);
    int_act = 1:timing(1);
    int_act = (int_act-1)*(activity(1)-mean(activity))/(timing(1)-1)+mean(activity);
    activities_int(1:timing(1),ind_n) = int_act;
    for ind_t=1:numel(activity)-1
        int_act = timing(ind_t):timing(ind_t+1);
        int_act = (int_act-timing(ind_t))*(activity(ind_t+1)-activity(ind_t))/(timing(ind_t+1)-timing(ind_t))+activity(ind_t);
        activities_int(timing(ind_t):timing(ind_t+1),ind_n) = int_act;
    end
    int_act = timing(end):duration;
    if numel(int_act)>1
        int_act = (int_act-timing(end))*(mean(activity)-activity(end))/(duration-timing(end))+activity(end);
    else
        int_act = activity(end);
    end
    activities_int(timing(end):duration,ind_n) = int_act;
end
end

function [signal_to_noise_ratio,maximo,lower_standard_deviation] = SNRChen_calculation(trace, percentile, max_threshold_factor, stim_time, panel,timing)
%compute SNR as in Chen et al. 2013 if there are stimuli or as in the GUI
%if not.
%max_threshold_factor controls the threshold for the trace points to be
%considered maxima. The larger max_threshold_factor, the higher the
%threshold. In the limit, the only possible maximum is the global one.
if stim_time==0
    lower_trace_th_aux = prctile(trace,percentile);
    lower_trace_index = trace<=lower_trace_th_aux;
    lower_trace_activity = trace(lower_trace_index);
    if std(lower_trace_activity)==0
        lower_standard_deviation = inf;
    else
        lower_standard_deviation = std(lower_trace_activity);
    end
    
    index1 = trace>max(trace)-(max(trace)-mean(lower_trace_activity))/max_threshold_factor;
    index2 = [false; diff(trace)>0];
    index3 = [diff(trace)<0; false];
    index = index1 & index2 & index3;
    if nnz(index)>0
        maximo = mean(trace(index));
    else
        maximo = max(trace);
    end
else
    %index for the baseline period
    lower_trace_index = 1:stim_time;
    lower_trace_activity = trace(lower_trace_index);
    %compute the (F-F0)/F0
    %trace = (trace-mean(lower_trace_activity))/(mean(lower_trace_activity)+eps);
    trace = 100*(trace-mean(lower_trace_activity))/(mean(lower_trace_activity)+eps);
    %get std on deltaF/f0 trace
    lower_trace_activity = trace(lower_trace_index);
    lower_standard_deviation = std(lower_trace_activity);
    trace_response = trace(stim_time+1:end);
    index1 = trace_response>max(trace_response)-(max(trace_response)-mean(lower_trace_activity))/max_threshold_factor;
    index2 = [true; diff(trace_response)>0];
    index3 = [diff(trace_response)<0; true];
    index = index1 & index2 & index3;
    if nnz(index)>0
        maximo = mean(trace_response(index));
    else
        maximo = max(trace_response);
    end
end

maximo = maximo-mean(lower_trace_activity);
signal_to_noise_ratio = maximo/(lower_standard_deviation+eps);
if isnan(signal_to_noise_ratio) || isinf(signal_to_noise_ratio)
    keyboard
end
if panel~=0
    timing = timing/1000;
    frame_rate = (1/(timing(2)-timing(1)));
    subplot(4,4,panel)
    hold on
    plot(timing,trace)
    plot(timing(find(index)+stim_time),trace(find(index)+stim_time),'+m','markerSize',10)
    plot(timing(lower_trace_index),lower_trace_activity,'+c','markerSize',10)
    if panel==1
        xlabel(gca,'time (s)')
        ylabel(gca,'fluorescence (%)')
        %legend(gca,'smart-line scan','raster scan','ground truth')
    end
    title([num2str(frame_rate,4) 'Hz | SNRChen:' num2str(signal_to_noise_ratio,4)])% ' DR: ' num2str(maximo,2) ' std: ' num2str(lower_standard_deviation,2)])
    set(gca,'fontSize',16)
    if stim_time~=0
        plot(timing(stim_time)*ones(1,2),[min(trace) max(trace)],'--','color',[.2 .2 .2])
    end
end
end

function [SNRMarco_scan, SNRMarco_max, SNRMarco_std] = SNRMarco_calculation(trace, spiketimes, scan_time, baseline_window, peak_window)
%computes SNR from transients associated to isolated spikes (spiketimes)

%baseline_window = baseline_window-5000;
SNRMarco_scan = zeros(1,numel(spiketimes));
SNRMarco_max = zeros(1,numel(spiketimes));
SNRMarco_std = zeros(1,numel(spiketimes));
for ind_spk=1:numel(spiketimes)
    %get DeltaF/F0
    baseline_ind = scan_time>=spiketimes(ind_spk)-baseline_window & scan_time<=spiketimes(ind_spk);
    baseline_trace = trace(baseline_ind);
    trace_deltaF_F0 = 100*(trace-mean(baseline_trace))/(mean(baseline_trace)+eps);
    
    %get baseline/transient values
    baseline_trace = trace_deltaF_F0(baseline_ind);
    transient_ind = scan_time>spiketimes(ind_spk) & scan_time<=spiketimes(ind_spk)+peak_window;
    transient_trace = trace_deltaF_F0(transient_ind);
    SNRMarco_scan(ind_spk) = max(transient_trace)/std(baseline_trace);
    SNRMarco_max(ind_spk) = max(transient_trace);
    SNRMarco_std(ind_spk) = std(baseline_trace);
    
    %plot trace
    %{
        figure
        hold on
        plot(scan_time,trace_deltaF_F0)
        plot(scan_time(baseline_ind),baseline_trace)
        plot(scan_time(transient_ind),transient_trace)
        plot([spiketimes(ind_spk) spiketimes(ind_spk)+peak_window],max(transient_trace)*ones(1,2),'--k')
        title([num2str(SNRMarco_scan(ind_spk)) ' | ' num2str(SNRMarco_max(ind_spk)) ' | ' num2str(SNRMarco_std(ind_spk))] )
        close all
    %}
end
SNRMarco_scan = mean(SNRMarco_scan);
SNRMarco_max = mean(SNRMarco_max);
SNRMarco_std = mean(SNRMarco_std);

end

function [axes_rasters,axes_stats] = figure_panels()
%panels for summary figure
left_margin = .075;
top_margin = .05;
small_panels_size = [.2 .35];
big_panels_size = [.3 .35];


axes_rasters(1) = axes('Position',[left_margin 1-top_margin-small_panels_size(2) small_panels_size(1) small_panels_size(2)]);
axes_rasters(2) = axes('Position',[left_margin 1-top_margin-small_panels_size(2)-(small_panels_size(2)+small_panels_size(2)/2) small_panels_size(1) small_panels_size(2)]);


axes_rasters(3) = axes('Position',[left_margin+small_panels_size(1)+small_panels_size(1)/2 1-top_margin-small_panels_size(2) small_panels_size(1) small_panels_size(2)]);
axes_rasters(4) = axes('Position',[left_margin+small_panels_size(1)+small_panels_size(1)/2 1-top_margin-small_panels_size(2)-(small_panels_size(2)+small_panels_size(2)/2) small_panels_size(1) small_panels_size(2)]);


% top_margin = top_margin + .025;
axes_stats(1) = axes('Position',[left_margin+2*small_panels_size(1)+small_panels_size(1)/2+big_panels_size(1)/4 ...
    1-top_margin-big_panels_size(2) big_panels_size(1) big_panels_size(2)]);
axes_stats(2) = axes('Position',[left_margin+2*small_panels_size(1)+small_panels_size(1)/2+big_panels_size(1)/4 ...
    1-top_margin-2*big_panels_size(2)-small_panels_size(2)/2 big_panels_size(1) big_panels_size(2)]);
end



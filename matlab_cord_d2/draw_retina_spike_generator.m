clear

%{
% --- draw using retina_spike_detector
filename = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/2Hz_vertical_2016Dec06/spike_Retina_Retina.mat';
a = load(filename);

a.times = a.times - min(a.times);
a.senders = a.senders - 1; % retina index are 2~1601 -> 1~1600 in matlab

num_step = max(a.times);
num_neuron = 1600;
data = zeros(num_neuron, num_step);

% Failed detectors
failed_detectors = load('/home/kfujii2/newNEST2/ht_model_pablo_based/data/2Hz_vertical_2016Dec06/retina_failed_detectors.mat');
failed_detectors.detectors = failed_detectors.detectors - 1;

for t=1:1:num_step
    pos = find(a.times==t);
    fired_idx = a.senders(pos);
    data(fired_idx,t) = 1;
    %data(failed_detectors.detectors,t) = 0.5;
end
%}

% --- draw using spike_times for spike_generators
%filename = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/2Hz_vertical_2016Dec07/retina_spike_times.mat';
%filename = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/transform_mp4/retina_spike_times.mat';
filename = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human/retina_spike_times.mat';
a = load(filename);
data = a.spike_times;
data = data';
num_step = size(data,2);

for i=1:1:num_step
    imagesc((reshape(data(:,i),40,40))');
    axis equal
    
    title( strcat('step ',num2str(i)) );
    pause(0.5)
end
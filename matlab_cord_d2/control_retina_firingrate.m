clear

%exp_dir = '2Hz_vertical_rate500';
%exp_dir = 'walking_human_frateTest';
exp_dir = 'moving_rectangle2';
%filename = strcat('/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/',...
filename = strcat('/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/',...
                  exp_dir,...  
                  '/spike_Retina_Retina.mat');
src_data = load(filename);

src_data.times = src_data.times - min(src_data.times) + 1;

% Put spike data into a matrix
num_neuron = max(src_data.senders) - min(src_data.senders) + 1;
duration = max(src_data.times) - min(src_data.times);
spike_mat = zeros(num_neuron,duration);

for t = 1:1:duration
    tmp_t_idx = find(src_data.times==t);
    tmp_n_idx = src_data.senders(tmp_t_idx);
    spike_mat(tmp_n_idx,t) = 1;
end

mean_val = mean(sum(spike_mat,1))
std_val = std(sum(spike_mat,1))
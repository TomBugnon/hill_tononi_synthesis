clear
%close all

%% Set parameters
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/2Hz_vertical_rate100/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/random_rate100/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human_rate100_dst/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human_rate100_scrbl_xy/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human_rate100_scrbl_t/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/moving_rectangle2_rate100_dst/';
dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/scramble_test_moving_rectangle2_rate20_dst/';
filename_v = strcat(dirname, 'spike_Vp_v_L4_exc.mat');
filename_h = strcat(dirname, 'spike_Vp_h_L4_exc.mat');

time_threshold_min = 500;
time_threshold_max = 5000;

time_window_min = 1; %[ms]
time_window_max = 1000; %[ms]
window_step = 1;


%% Load data
% Output from NEST uses sparse-expression, so I converted to matrix
% spike_mat
% --- row: neuron index
% --- col: time sequence
% --- If i-th neuron fired at time j, spike_mat(i,j)==1, otherwise 0
% --- Note that time is limited from "time_window_min" to "time_window_max"
% --- so, j-th column in spike_mat represents the time of (j + time_window_min) in NEST simulator.
spike_mat_v = convert_nestspike2mat(filename_v,time_threshold_min,time_threshold_max);
spike_mat_h = convert_nestspike2mat(filename_h,time_threshold_min,time_threshold_max);

spike_mat = [spike_mat_v; spike_mat_h];
%spike_mat = spike_mat_v;

%% Preparation
% Ignore neurons which did not fired at all during experiments
ignore_idx = find( sum(spike_mat,2)==0 );
spike_mat(ignore_idx,:) = [];

% Define some parameters and containers
%duration = time_threshold_max - time_threshold_min; % length of data
%num_neuron = size(spike_mat,1); % the number of neurons
duration = time_threshold_max - time_threshold_min; % length of data
num_neuron = size(spike_mat,1); % the number of neurons

% --- The following valuables are containers to keep data at every step in the main loop.
% --- Because the number of steps depends on "time_window", 
% --- I didn't declare the size of those containers in advance.
d2_all = []; % a container to keep D2 value
threshold_all = []; % a container to keep threshold of firing rate to separate states
mean_p = []; % mean value of probability being 0/1 calculated from probability of all neurons
std_p = []; % std value of probability being 0/1 calculated from probability of all neurons
firing_hist_all = [];
median_each_all=[];

x_val = []; % xlabel value
x_time = [];
x_str = {}; % xlabel string
loop = 0; % someting to make x_val and x_str

%% Main loop
% In order to find the time resolution which gives maximumn state differentiation(D2),
% --- 1. Calculate firing rate within a time window
% --- 2. Find thresold value(firing_rate) to separate states(single neuron) into 2 states.
% --- 3. Calculate state probability for each neuron using the threshold value
% --- 4. Calculate state differentiation(d2)


time_steps = time_window_min:window_step:time_window_max;
thresholds = 1:1:100;
%thresholds = 1:2:100;

all_d2 = zeros(length(thresholds), length(time_steps));
for firing_index = 1:length(thresholds)
    
    firing_threshold = thresholds(firing_index);
    
    fprintf('Computing D2s for threshold == %d\n', firing_threshold);
    for time_index = 1:length(time_steps)
        
        time_window = time_steps(time_index);
        
        % 1. Calulate firing rate within a time-window for each neuron
        num_step = floor(duration/time_window);
        firing_rate = zeros(size(spike_mat,1),num_step);
        s = 1;
        e = 1;
        for step = 1:1:num_step
            e = s + time_window - 1;
            time_bin = spike_mat(:,s:e);
            %firing_rate(:,step) = sum(time_bin,2)/time_window; % #_of_spikes / time_window[ms]
            firing_rate(:,step) = sum(time_bin,2)/(time_window*0.001); % #_of_spikes / time_window[sec]
            s = e + 1;
        end

        % 3. Calculate state probability for each neuron using the threshold value
        % --- Segregate states using threshold
        states = zeros(num_neuron,num_step);
        idx = find(firing_rate>firing_threshold);
        states(idx) = 1;

        % --- Calculate state probability
        % --- if states(i)==0 then count as state1 -> probability p(i,1)
        % --- if states(i)==1 then count as state2 -> probability p(i,2)
        p = zeros(num_neuron,2);
        p(:,2) = sum(states,2) ./ num_step; % more active state
        p(:,1) = (repmat(num_step,num_neuron,1)-sum(states,2)) / num_step; % less active state

        % 4. Calculate state differentiation(d2)
        % --- p_positive only contains probabilties such that both of state 0/1 have positive values
        % --- the reason I ignore other parobabilities is,
        % --- if probability for state 0(or1) is 1.0, the entropy will be zero,
        % --- but the calculation will result in nan ( 0*log2(0)=nan )
        p_positive = p;
        state1_zero_idx = find(p_positive(:,1)==0);
        state2_zero_idx = find(p_positive(:,2)==0);
        remove_idx = [state1_zero_idx; state2_zero_idx];
        p_positive(remove_idx,:) = [];

    %     entropy_each = -1 * ( p_positive(:,1).*log2(p_positive(:,1)) +...
    %                           p_positive(:,2).*log2(p_positive(:,2)) );
    %     d2 = sum(entropy_each);

        variance_each = p_positive(:,1).*(1-p_positive(:,1)) + p_positive(:,2).*(1-p_positive(:,2));
        d2 = sum(variance_each);

        d2 = d2/num_neuron;

        % For analysis 
        % Keep the parameters which gives the maximumn D2
        if d2 > max(d2_all)
            max_d2 = d2;
            max_p_positive = p_positive;
            max_p = p;
            max_resolution = time_window;
            max_step = loop;
            max_threshold = firing_threshold;
            %max_median_each = median_each;
            max_firing_rate = firing_rate;
        end

    %     if any(ismember(time_window, [53,55,267]))
    %         figure
    %         histogram(p(:,1),20,'FaceColor','k')
    %         hold on
    %         histogram(p(:,2),20,'FaceColor','r')
    %         title(time_window)
    %     end

        % Keep other information
        % --- Note that the following codes make processing slower.
        % --- If you need faster processing, deleting this section might be a solution.
        d2_all = [d2_all; d2];
        threshold_all = [threshold_all; firing_threshold];
        mean_p = [mean_p; mean(p,1)];
        std_p = [std_p; std(p,1)];
        h = hist(reshape(firing_rate,1,numel(firing_rate)), [0:1:200]);
        firing_hist_all = [firing_hist_all; h];
        %median_each_all = [median_each_all; median_each'];

        % For a graph 
        loop = loop + 1;
        x_val = [x_val; loop];
        x_time = [x_time; time_window];
        x_str = [x_str; num2str(time_window)];

        % keyboard
        all_d2(firing_index, time_index) = d2;
    end
end

%%
figure
imagesc(time_steps,thresholds,all_d2);
title_str = strsplit(dirname, '/');
title(strrep(title_str{end-1}, '_', '\_'))
ylabel('Threshold (Hz)');
xlabel('Window (ms)');
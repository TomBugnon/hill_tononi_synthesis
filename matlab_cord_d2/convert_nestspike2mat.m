function spike_mat = convert_nestspike2mat(filename, time_threshold_min, time_threshold_max)

    src_data = load(filename);

    % Remove unnecessory part of the data
    early_idx = find(src_data.times<time_threshold_min);
    src_data.senders(early_idx) = [];
    src_data.times(early_idx) = [];

    late_idx = find(src_data.times>time_threshold_max);
    src_data.senders(late_idx) = [];
    src_data.times(late_idx) = [];

    % Convert index from NEST to MATLAB
    src_data.senders = src_data.senders - min(src_data.senders) + 1;
    src_data.times = src_data.times - min(src_data.times) + 1;

    % Put spike data into a matrix
    num_neuron = max(src_data.senders);
    duration = time_threshold_max - time_threshold_min;
    %duration = max(src_data.times) - min(src_data.times);
    spike_mat = zeros(num_neuron,duration);
    for t = 1:1:duration
        tmp_t_idx = find(src_data.times==t);
        tmp_n_idx = src_data.senders(tmp_t_idx);
        spike_mat(tmp_n_idx,t) = 1;
    end
    
end
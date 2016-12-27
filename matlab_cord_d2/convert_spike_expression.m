function spike = convert_spike_expression( sparse_mat, num_draw_neurons, resolution )

    % Input format
    % --- # of rows: neuron index
    % --- # of cols: different for each neuron depening on how many times that neuron fired
    % --- Data are float value which represent spike timing
    % --- Note: 'sparse_mat' is cell struct
    %
    % Output format 
    % --- # of rows: neuron index
    % --- # of cols: time duration
    % --- Data are 0/1 binary value which represent whether a neuron fired or not at a time point
    % --- Note: 'sparse_mat' is double mat

    % Arguments
    % 1. sparse_mat (See "Input format" above.)
    % 2. num_draw_neurons
    % --- The number of neurons which project to V1
    % --- This is not necessarily the same as num_all_neurons(=9,000) 
    % --- Specifically, num_draw_neurons = length( "gid in LGN_information.csv" );
    % 3. resolution
    % --- Allen data has micro second time resolution.
    % --- resolution = 1 -> Time resolution after processing will be 1 sec
    % --- resolution = 10 -> Time resolution after processing will be 0.1 sec


    % Time length
    max_time = ceil( max( cell2mat(sparse_mat)*resolution ) );
    min_time = ceil( min( cell2mat(sparse_mat)*resolution ) );
    duration = max_time - min_time + 1;

    spike = zeros(num_draw_neurons,duration);

    % Set spike timing data to 'spike'
    for n=1:1:num_draw_neurons

        % Spike timing is float value in Allen, so round it.
        % Start timing is not necessarily 1.0, so subtract 'min_time'
        spike_time = ceil( sparse_mat{n}*resolution );
        spike_time = spike_time - min_time + 1;

        % Take unique of 'spike_time'
        % 'spike_time' can be redundant because of ceil()
        spike_time = unique(spike_time);
        spike(n,spike_time) = 1;

    end
    
end
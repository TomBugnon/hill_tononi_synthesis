function recorded_activity = convert_nestrecorder2mat(filename, time_threshold_min, time_threshold_max, fields)
    % fields is cellarray containing the fields of interest.
    % Returns structure with neuron*time array for each field
    
    recorded_activity = struct();
    
    
    src_data = load(filename);
    
    
    
    
    
    
    
    % Remove unnecessory part of the data
    early_idx = find(src_data.times<time_threshold_min);
    srcfields = fieldnames(src_data);
    for fieldnum = 1:numel(srcfields)
        src_data.(srcfields{fieldnum})(early_idx) = [];
        
    end
    
    
    
    if time_threshold_max < 0
        % round the maximum by 100 ms, assume the network rund in multiples
        % of 100 ms
        time_threshold_max = round(max(src_data.times)/100)*100;
    else
        late_idx = find(src_data.times>time_threshold_max);
        srcfields = fieldnames(src_data)
        for fieldnum = 1:numel(srcfields)
            src_data.(srcfields{fieldnum})(late_idx) = [];
            
        end
    end
    
    
    
    % Convert index from NEST to MATLAB
    src_data.senders = src_data.senders - min(src_data.senders) + 1;
    src_data.times = src_data.times - min(src_data.times) + 1;
    
    
    %Prepare output fields
%     num_neuron = max(src_data.senders)- min(src_data.senders);
    num_neuron = max(src_data.senders);

    duration = time_threshold_max - time_threshold_min;
    for fieldnumb = 1:numel(fields)
        recorded_activity.(fields{fieldnumb}) = zeros(num_neuron,duration);
    end
    
    %duration = max(src_data.times) - min(src_data.times);
    
    %For each field, put data in array that is in field of main struct
    for t = 1:1:duration
        tmp_t_idx = find(src_data.times==t);
        tmp_n_idx = src_data.senders(tmp_t_idx);
        for fieldnumb = 1:numel(fields)
            recorded_activity.(fields{fieldnumb})(tmp_n_idx,t) = src_data.(fields{fieldnumb})(tmp_t_idx);
        end
    end
end
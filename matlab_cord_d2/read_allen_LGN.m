function tmp_trial_data = read_allen_LGN (filename, trial_idx, num_trials)
    
    % Arguments
    % 1. filename
    % --- filename which has LGN spike data.
    % --- This file has [num_all_neurons(=9,000) * num_trials] lines
    % --- If num_trials is 10, rows 1-10 represents spike_timing of LGN_neuron 0 in trial 1-10
    % 2. trial_idx
    % --- Index of trial of interest.
    % --- This should be int value (Range is 1 ~ num_trials)
    % 3. num_trials
    % --- 20 for Movie data,
    % --- 10 for grating data.
    % --- Check the number of rows of the file. Ex: wc -l filename
    % --- num_trials = (the number of rows) / (num_all_neurons(=9000))
    
    fid = fopen(filename);
    disp( strcat('Reading: ', filename) )

    % Prepare container 
    % Use cell matrix without difining size because
    % 1. the number of rows is not clear
    % 2. the number of colomns is not constant
    all_data = {};
    num_lines = 1;

    % Get the first line of the file
    tmp_line = fgetl(fid);

    % Read each line and append data to 'all_data'
    while ischar(tmp_line)

        % To see the progress of this code.
        if mod(num_lines,2000)==0
            num_lines
        end

        % if the first character is not a number, remove it.
        if tmp_line(1)==' '
            tmp_line(1)='';
        end

        % Split tmp_line to several parts
        % Delimiter is space(' ') in Allen's LGN data
        tmp_string = strsplit(tmp_line,' ');

        % Convert string to number
        tmp_data = zeros(size(tmp_string));
        for i=1:1:length(tmp_string)
            tmp_data(i) = str2double(tmp_string{i}); 
        end

        % Append 'tmp_data' as the last element of 'all_data'
        all_data{num_lines} = tmp_data;

        % Prepare for the next loop
        % When the file ends 'tmp_line' will be -1.
        % That results in exiting from while-loop (because ischar(tmp_line)=False)
        num_lines = num_lines + 1;
        tmp_line = fgetl(fid);

    end

    fclose(fid)
    
    num_all_neurons = 9000;
    num_lines = num_all_neurons*num_trials;
    tmp_trial_lines = trial_idx:num_trials:num_lines;
    tmp_trial_data = all_data(tmp_trial_lines);
    
end
function scrambled_connection = scramble_post_neuron(original_data,post_min,post_max)

    tmp_ori = original_data(ismember(original_data(:,2),[post_min:1:post_max]), :);    
    new_post = randi([L4_h_exc_min, L4_h_exc_max], size(tmp_ori,1), 1);
    scrambled_connection = tmp_ori;
    scrambled_connection(:,2) = new_post; % replace post neuron index
end
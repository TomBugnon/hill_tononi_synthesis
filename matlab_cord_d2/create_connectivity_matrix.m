function cmat = create_connectivity_matrix(filename, target_min, target_max)

    %dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/2Hz_vertical_rate100_alldata/';
    %filename = 'connection_Vp_v_AMPA_syn.dat';
    %src_data = load( strcat(dirname,filename) );
    src_data = load(filename);

    %target_min = 20806; % minimum index of V1 vertical population
    %target_max = 35205; % maximum index of V1 vertical population
    num_neurons = target_max - target_min + 1;
    cmat = zeros(num_neurons,num_neurons); % connectivity matrix

    tmp_pre_rows = ismember( src_data(:,1), [target_min:1:target_max] );
    tmp_post_rows = ismember( src_data(:,2), [target_min:1:target_max] );
    rows = find( (tmp_pre_rows+tmp_post_rows)==2 );

    pre_idx = src_data(rows,1) - target_min + 1;
    post_idx = src_data(rows,2) - target_min + 1;

    connected_idx = pre_idx + num_neurons*(post_idx-1);

    cmat(connected_idx) = 1;
end
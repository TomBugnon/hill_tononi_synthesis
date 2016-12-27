% stat_NESTconnection()
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/2Hz_vertical_rate100_alldata/';
dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/scramble_test2_moving_rectangle2_rate20_dst/';
%filename_GABAA = strcat(dirname,'connection_Tp_GABA_A_syn.dat');
%all_data = load(filename_GABAA);
filename_AMPA = strcat(dirname,'connection_Tp_AMPA_syn.dat');
%filename_AMPA = strcat(dirname,'connection_Retina_AMPA_syn.dat');
all_data = load(filename_AMPA);

% format of all_data
% --- column1: pre neuron gid
% --- column2: post neuron gid
% --- column3: weight
% --- column4: delay
% --- column5: dx
% --- column6: dy
% --- column7: dz(if the layer has 3D structure)

layer_size = 40*40;

%% Set neuron index for each populaiton
% --- Horizontal L23,L4,L56 / exc,inh
horizontal_head = 6405;
L23_h_exc_min  = horizontal_head;
L23_h_exc_max  = L23_h_exc_min + (layer_size*2) - 1;
L23_h_inh_min  = L23_h_exc_max + 1;
L23_h_inh_max  = L23_h_inh_min + layer_size - 1;
%
L4_h_exc_min  = L23_h_inh_max + 1;
L4_h_exc_max  = L4_h_exc_min + (layer_size*2) - 1;
L4_h_inh_min  = L4_h_exc_max + 1;
L4_h_inh_max  = L4_h_inh_min + layer_size - 1;
%
L56_h_exc_min = L4_h_inh_max + 1;
L56_h_exc_max = L56_h_exc_min + (layer_size*2) - 1;
L56_h_inh_min = L56_h_exc_max + 1;
L56_h_inh_max = L56_h_inh_min + layer_size - 1;
%
% --- Vertical L23,L4,L56 / exc,inh
vertical_head = 20806;
L23_v_exc_min = vertical_head;
L23_v_exc_max = L23_v_exc_min + (layer_size*2) - 1;
L23_v_inh_min = L23_v_exc_max + 1;
L23_v_inh_max = L23_v_inh_min + layer_size - 1;
%
L4_v_exc_min  = L23_v_inh_max + 1;
L4_v_exc_max  = L4_v_exc_min + (layer_size*2) - 1;
L4_v_inh_min  = L4_v_exc_max + 1;
L4_v_inh_max  = L4_v_inh_min + layer_size - 1;
%
L56_v_exc_min = L4_v_inh_max + 1;
L56_v_exc_max = L56_v_exc_min + (layer_size*2) - 1;
L56_v_inh_min = L56_v_exc_max + 1;
L56_v_inh_max = L56_v_inh_min + layer_size - 1;

%% Count the number of connections from Tp to each layers
% Horizontal
num_L23_h_exc = sum( ismember(all_data(:,2),[L23_h_exc_min:1:L23_h_exc_max]) );
num_L23_h_inh = sum( ismember(all_data(:,2),[L23_h_inh_min:1:L23_h_inh_max]) );
num_L4_h_exc  = sum( ismember(all_data(:,2),[L4_h_exc_min:1:L4_h_exc_max]) );
num_L4_h_inh  = sum( ismember(all_data(:,2),[L4_h_inh_min:1:L4_h_inh_max]) );
num_L56_h_exc = sum( ismember(all_data(:,2),[L56_h_exc_min:1:L56_h_exc_max]) );
num_L56_h_inh = sum( ismember(all_data(:,2),[L56_h_inh_min:1:L56_h_inh_max]) );

% vertical
num_L23_v_exc = sum( ismember(all_data(:,2),[L23_v_exc_min:1:L23_v_exc_max]) );
num_L23_v_inh = sum( ismember(all_data(:,2),[L23_v_inh_min:1:L23_v_inh_max]) );
num_L4_v_exc  = sum( ismember(all_data(:,2),[L4_v_exc_min:1:L4_v_exc_max]) );
num_L4_v_inh  = sum( ismember(all_data(:,2),[L4_v_inh_min:1:L4_v_inh_max]) );
num_L56_v_exc = sum( ismember(all_data(:,2),[L56_v_exc_min:1:L56_v_exc_max]) );
num_L56_v_inh = sum( ismember(all_data(:,2),[L56_v_inh_min:1:L56_v_inh_max]) );

%%
% --- Horizontal
% L4 exc
ori_data = all_data(ismember(all_data(:,2),[L4_h_exc_min:1:L4_h_exc_max]), 1:4);    
new_post = randi([L4_h_exc_min, L4_h_exc_max], size(ori_data,1), 1);
Tp_L4_h_exc = ori_data;
Tp_L4_h_exc(:,2) = new_post; % replace post neuron index
% L4 inh
ori_data = all_data(ismember(all_data(:,2),[L4_h_inh_min:1:L4_h_inh_max]), 1:4);    
new_post = randi([L4_h_inh_min, L4_h_inh_max], size(ori_data,1), 1);
Tp_L4_h_inh = ori_data;
Tp_L4_h_inh(:,2) = new_post; % replace post neuron index

% L56 exc
ori_data = all_data(ismember(all_data(:,2),[L56_h_exc_min:1:L56_h_exc_max]), 1:4);
new_post = randi([L56_h_exc_min, L56_h_exc_max], size(ori_data,1), 1);
Tp_L56_h_exc = ori_data;
Tp_L56_h_exc(:,2) = new_post;
% L56 inh
ori_data = all_data(ismember(all_data(:,2),[L56_h_inh_min:1:L56_h_inh_max]), 1:4);
new_post = randi([L56_h_inh_min, L56_h_inh_max], size(ori_data,1), 1);
Tp_L56_h_inh = ori_data;
Tp_L56_h_inh(:,2) = new_post;

% --- Vertical exc
% L4 exc
ori_data = all_data(ismember(all_data(:,2),[L4_v_exc_min:1:L4_v_exc_max]), 1:4);    
new_post = randi([L4_v_exc_min, L4_v_exc_max], size(ori_data,1), 1);
Tp_L4_v_exc = ori_data;
Tp_L4_v_exc(:,2) = new_post; 
% L4 inh
ori_data = all_data(ismember(all_data(:,2),[L4_v_inh_min:1:L4_v_inh_max]), 1:4);    
new_post = randi([L4_v_inh_min, L4_v_inh_max], size(ori_data,1), 1);
Tp_L4_v_inh = ori_data;
Tp_L4_v_inh(:,2) = new_post; 

% L56 exc
ori_data = all_data(ismember(all_data(:,2),[L56_v_exc_min:1:L56_v_exc_max]), 1:4);
new_post = randi([L56_v_exc_min, L56_v_exc_max], size(ori_data,1), 1);
Tp_L56_v_exc = ori_data;
Tp_L56_v_exc(:,2) = new_post;
% L56 inh
ori_data = all_data(ismember(all_data(:,2),[L56_v_inh_min:1:L56_v_inh_max]), 1:4);
new_post = randi([L56_v_inh_min, L56_v_inh_max], size(ori_data,1), 1);
Tp_L56_v_inh = ori_data;
Tp_L56_v_inh(:,2) = new_post;

scrambled_connections = [Tp_L4_h_exc;...
                         Tp_L4_h_inh;...
                         Tp_L56_h_exc;...
                         Tp_L56_h_inh;...
                         Tp_L4_v_exc;...
                         Tp_L4_v_inh;...
                         Tp_L56_v_exc;...
                         Tp_L56_v_inh];
dlmwrite('scrambled_connection_Tp_Cortex.dat', 'delimiter', ' ');
save('scrambled_connection_Tp_Cortex.mat', 'scrambled_connections');

%% Calc connection probabilities
% % Horizontal
% p_L23_h_exc = num_L23_h_exc / (layer_size*2);
% p_L23_h_inh = num_L23_h_inh / layer_size;
% p_L4_h_exc  = num_L4_h_exc  / (layer_size*2);
% p_L4_h_inh  = num_L4_h_inh  / layer_size;
% p_L56_h_exc = num_L56_h_exc / (layer_size*2);
% p_L56_h_inh = num_L56_h_inh / layer_size;
% 
% % vertical
% p_L23_v_exc = num_L23_v_exc / (layer_size*2);
% p_L23_v_inh = num_L23_v_inh / layer_size;
% p_L4_v_exc  = num_L4_v_exc  / (layer_size*2);
% p_L4_v_inh  = num_L4_v_inh  / layer_size;
% p_L56_v_exc = num_L56_v_exc / (layer_size*2);
% p_L56_v_inh = num_L56_v_inh / layer_size;
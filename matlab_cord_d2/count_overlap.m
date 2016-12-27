clear

%% Set parameters
dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/2Hz_vertical_rate100_alldata/';
filename = 'connection_Vp_v_AMPA_syn.dat';

V1_v_all_min = 20806; % minimum index of V1 vertical population
V1_v_all_max = 35205; % maximum index of V1 vertical population
num_exc = (40*40)*2; % the number of excitatory neurons in a single layer (L23 for example)
num_inh = (40*40); % the number of inhibitory neurons in a single layer (L23 for example)

%% Create connectivity matrix
    
L4_v_exc_min = V1_v_all_min + num_exc + num_inh; % minimum index of V1-Layer4-excitatory neurons
%L4_v_exc_max = L4_v_exc_min + num_exc; % maximum index of V1-Layer4-excitatory neurons
L4_v_exc_max = L4_v_exc_min + (40*40); % maximum index of V1-Layer4-excitatory neurons

cmat = create_connectivity_matrix( strcat(dirname,filename), L4_v_exc_min, L4_v_exc_max );

%% Count overlaps

%pre_overlap = sum(cmat,2);
%post_overlap = sum(cmat,1);

if size(cmat,1)~=size(cmat,2)
    disp('check size of cmat. # of rows and # of col should be equal.')
else
    num_neuron = size(cmat,1);
end
mechanism_pre = cell(num_neuron,1);
mechanism_post = cell(num_neuron,1);

tic
for i=1:1:num_neuron
    
    %post_idx = i;
    tmp_post_mech = find(cmat(:,i)>0);
    if length(tmp_post_mech)>2
        mechanism_post{i} = tmp_post_mech;
    end
    
    %pre_idx = i;
    tmp_pre_mech = find(cmat(i,:)>0);
    if length(tmp_pre_mech)>2
        mechanism_pre{i} = tmp_pre_mech;
    end    
    
end
toc

i=1
target_mech = mechanism_pre{1};



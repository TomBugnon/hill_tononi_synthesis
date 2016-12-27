%% Parameters

%--- Movie
num_trials = 20;
trial_idx = 1; % 0 ~ num_trials
%struct_filename = 'data/allen_tmp/spk_lgn/mov_2_LGN_spk.dat';
random_filename = 'data/allen_tmp/spk_lgn/mov_2_scrbl_t_LGN_spk.dat';

%--- grating
%num_trials = 10;
%trial_idx = 1; % 0 ~ num_trials
struct_filename = 'data/allen_tmp/spk_lgn/grating_7_LGN_spk.dat';

%% Load LGN information
% This file is made by remving the following elements from 'LGN_cells_information.csv'
% - The first row (explanation of each column data)
% - LGN neuron type column (transient ON/OFF/ON-OFF)
% - Synaptic type column (LGN exc)
LGN_info = load('data/allen_tmp/spk_lgn/LGN_cells_information_copy.csv');

% Select data only needed in this code
% [Remove] 1st column: target V1 neuron index
% [Select] 2nd column: source LGN neuron index -> 'gid'
% [Select] 3rd column: x position of a LGN neuron -> 'x'
% [Select] 4rd column: y position of a LGN neuron -> 'y'
% [Remove] 5th column: # of synapses from a source LGN neuron to a target V1 neuron
LGN_info = LGN_info(:,2:4);
% Data is redundant, so take unique rows of LGN_info
LGN_info = unique(LGN_info,'rows');

gid = LGN_info(:,1);
gid = gid + 1; % convert index: In allen's format, index starts from one, but matlab index starts from 1.
x = LGN_info(:,2);
y = LGN_info(:,3);


%% Read LGN_*_spk.dat data
%struct_data_sparse = read_allen_LGN(struct_filename, trial_idx, num_trials);
struct_data_sparse = read_allen_LGN(struct_filename, trial_idx, 10);
random_data_sparse = read_allen_LGN(random_filename, trial_idx, num_trials);

%% Convert spike data 
num_draw_neurons = length(gid); 
resolution = 1;
struct_spike = convert_spike_expression(struct_data_sparse, num_draw_neurons, resolution);
random_spike = convert_spike_expression(random_data_sparse, num_draw_neurons, resolution);

%{
%% Draw 2D mapping
s_duration = size(struct_spike,2);
r_duration = size(random_spike,2);
duration = min(s_duration,r_duration);
figure;
for t=1:1:duration
    struct_fired_neurons = find(struct_spike(:,t)==1);
    random_fired_neurons = find(random_spike(:,t)==1);
    
    scatter(x(struct_fired_neurons), y(struct_fired_neurons), 200, '.', 'r');
    hold on
    scatter(x(random_fired_neurons), y(random_fired_neurons), 200, '.', 'b');
    axis equal
    grid on
    legend('struct', 'random')
    xlim([35 205])
    ylim([10 110])
    title( strcat('t = ',num2str(t)) )
    hold off
    pause(0.1)
end
%}

%% distance_analysis
if size(struct_spike,2)~=size(random_spike,2)
    disp('time length is not the same.')
    dif_length = abs(size(struct_spike,2) - size(random_spike,2));

    if size(struct_spike,2) > size(random_spike,2)
        struct_spike(:,1:dif_length) = [];
        disp(strcat('struct_spike is longer, removed ', num2str(dif_length), 'columns'))
    end
    if size(struct_spike,2) < size(random_spike,2)
        random_spike(:,1:dif_length) = [];
        disp(strcat('random_spike is longer, removed ', num2str(dif_length), 'columns'))
    end
end
time_length = size(struct_spike,2);
    
struct_num_all_spike = length(struct_spike==1);
random_num_all_spike = length(random_spike==1);
struct_distance = zeros(struct_num_all_spike,1);
random_distance = zeros(random_num_all_spike,1);

struct_num_firing = zeros(time_length,1);
random_num_firing = zeros(time_length,1);

struct_time_firing = zeros(struct_num_all_spike,1);
random_time_firing = zeros(random_num_all_spike,1);

ss = 1; % ss: Struct Start 
se = ss; % se: Struct End
rs = 1; % rs: Random Start 
re = rs; % se: Random End

for t=1:1:time_length
    
    % structured
    struct_isfired = (struct_spike(:,t)==1);
    struct_fired_x = x(struct_isfired);
    struct_fired_y = y(struct_isfired);
    struct_tmp_d = sqrt( struct_fired_x.^2 + struct_fired_y.^2 );
    
    se = ss + length(struct_tmp_d) - 1;
    struct_distance(ss:se) = struct_tmp_d;
    struct_num_firing(t) = length(struct_tmp_d);
    struct_time_firing(ss:se) = repmat(t, length(struct_tmp_d),1);
    ss = se + 1;
    
    % random
    random_isfired = (random_spike(:,t)==1);
    random_fired_x = x(random_isfired);
    random_fired_y = y(random_isfired);
    random_tmp_d = sqrt( random_fired_x.^2 + random_fired_y.^2 );

    re = rs + length(random_tmp_d) - 1;
    random_distance(rs:re) = random_tmp_d;
    random_num_firing(t) = length(random_tmp_d);
    random_time_firing(rs:re) = repmat(t, length(random_tmp_d),1);
    rs = re + 1;
    
end

figure

t_step = 2;
n_step = 2;

struct_spike_draw = struct_spike(1:n_step:end,:);
struct_spike_draw = struct_spike_draw(:,1:t_step:end);

random_spike_draw = random_spike(1:n_step:end,:);
random_spike_draw = random_spike_draw(:,1:t_step:end);

subplot(1,3,1)
h1 = pcolor(struct_spike_draw);
set(h1,'EdgeColor','none')
title(['Structured LGN'])
%xlim([0 1000])
subplot(1,3,2)
h2 = pcolor(random_spike_draw);
set(h2,'EdgeColor','none')
title(['Random LGN' ])
%xlim([0 1000])

subplot(1,3,3)
nbins = [0:50:time_length];
histogram(struct_time_firing,nbins) %keiko
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
%histogram(R) %original
histogram(random_time_firing,nbins) %keiko
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);

legend('structured','random')

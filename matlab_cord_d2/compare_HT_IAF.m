clear

hard_min_time = 500.0;
%hard_min_time = 0.0;
hard_max_time = 5000.0;

%hard_min_time = 0;
%hard_max_time = 0;

data_dir = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency';
% data_dir = '/home/kfujii2/newNEST2/ht_model_pablo_based/data';
random_data = '/random_2016Dec07_rate05';
structured_data = '/2Hz_vertical_2016Dec07_rate05';

% population_name = 'Retina'; plot_name = population_name;
%populations = {'Vp_v_L4_exc', 'Vp_h_L4_exc'}; plot_name = 'All L4pyr'; %original
%populations = {'Vp_v_L4_exc'}; plot_name = 'V L4pyr';
%populations = {'Vp_h_L4_exc'}; plot_name = 'H L4pyr';
%populations = {'Tp_Tp_exc'}; plot_name = 'Tp';

struct = load([data_dir structured_data '/spike_' populations{1} '.mat']);
struct.senders = struct.senders - min(struct.senders) + 1;
random = load([data_dir random_data '/spike_' populations{1} '.mat']);
random.senders = random.senders - min(random.senders) + 1;

%for population_name = populations(2:end) % original
for i = 2:1:size(populations,2) % keiko

    %--- original
    % file_to_load = ['/spike_' population_name{1}]; 
    %--- keiko
    population_name = populations{i};
    file_to_load = ['/spike_' population_name];

    s = load([data_dir structured_data file_to_load '.mat']);
    s.senders = s.senders - min(s.senders) + 1;

    struct.senders = [struct.senders  s.senders+max(struct.senders)+1];
    struct.times = [struct.times s.times];

    r = load([data_dir random_data file_to_load '.mat']);
    r.senders = r.senders - min(r.senders) + 1;
    
    random.senders = [random.senders r.senders+max(random.senders)+1];
    random.times = [random.times r.times];
    
end

% Convert neuon index
% Allen format starts from 0, matlab starts from 1
struct.senders = struct.senders + 1;
random.senders = random.senders + 1;


time_grain = 1; % mili seconds
% time_grain = 1000; %micro senconds

%convert time scale
struct.times = int32(time_grain * struct.times);
random.times = int32(time_grain * random.times);

% apply threshold
if hard_max_time ~= 0
    too_late = struct.times > hard_max_time;
    struct.times(too_late) = [];
    struct.senders(too_late) = [];
    
    too_late = random.times > hard_max_time;
    random.times(too_late) = [];
    random.senders(too_late) = [];
end

if hard_min_time ~= 0
    too_soon = struct.times < hard_min_time;
    struct.times(too_soon) = [];
    struct.senders(too_soon) = [];
    
    too_soon = random.times < hard_min_time;
    random.times(too_soon) = [];
    random.senders(too_soon) = [];
end


maxt = int32((100*time_grain)*(ceil(max(struct.times)/(100*time_grain))));
if maxt ~= int32((100*time_grain)*(ceil(max(random.times)/(100*time_grain))))
    disp('random and structured do not have the same duration!')
end

% discard first part of simulation
mint = int32(time_grain*min(struct.times));
if mint ~= int32(time_grain*min(random.times))
    disp('random and structured do not have the same startnig point!')
end

maxn = int32(100*(ceil(max(struct.senders)/100)));
if maxn ~= int32(100*(ceil(max(random.senders)/100)))
    disp('random and structured do not have the same # of neurons!')
end

minn = int32(min(struct.senders));
if minn ~= int32(min(random.senders))
    disp('random and structured do not have the same startnig point!')
end

struct.senders = int32(struct.senders);
random.senders = int32(random.senders);

fprintf('\n\nEstimated memory needed : %2.2f Gb\n\n', (2*maxt*maxn*8) / (2^30))


%% plot
%N_neuron = length(populations_idx);
N_neuron = max(random.senders);
N_sample = hard_max_time - hard_min_time + 1;
struct_2D_all = zeros(N_neuron, N_sample);
random_2D_all = zeros(N_neuron, N_sample);

column = 1;
for t = hard_min_time:1:hard_max_time
    tmp_struct_rows = find(struct.times==t);
    tmp_struct_senders = struct.senders(tmp_struct_rows);
    struct_2D_all(tmp_struct_senders,column) = 1;

    tmp_random_rows = find(random.times==t);
    tmp_random_senders = random.senders(tmp_random_rows);
    random_2D_all(tmp_random_senders,column) = 1;

    column = column + 1;
end

draw_neuron_idx = 1:15:N_neuron;
draw_time_idx = 1:2:N_sample;

struct_2D = struct_2D_all(draw_neuron_idx,draw_time_idx);
random_2D = random_2D_all(draw_neuron_idx,draw_time_idx);

%{
% downsample to visualitaion (screen resolution)

sizeN = 1600; %original
%sizeT = 2001; %original
%sizeN = 3200;
sizeT = 2001;

downsampleN = ceil(double(maxn-minn+1)/double(sizeN)); 
downsampleT = ceil(double(maxt-mint+1)/double(sizeT)); 

struct_2D = zeros(sizeN, sizeT);
random_2D = struct_2D;

for is = 1:length(struct.senders)
    struct_2D(ceil(double(struct.senders(is)-minn+1)/double(downsampleN)), ceil(double(struct.times(is)-mint+1)/double(downsampleT))) = 1; 
end

for is = 1:length(random.senders)
    random_2D(ceil(double(random.senders(is)-minn+1)/double(downsampleN)), ceil(double(random.times(is)-mint+1)/double(downsampleT))) = 1; 
end


% struct_2D = zeros(maxn-minn+1, maxt-mint);
% random_2D = struct_2D;
% 
% for is = 1:length(struct.senders)
%     struct_2D(struct.senders(is)-minn+1, struct.times(is)-mint+1) = 1; 
% end
% 
% for is = 1:length(random.senders)
%     random_2D(random.senders(is)-minn+1, random.times(is)-mint+1) = 1;
% end
%}

figure


subplot(1,3,1)
h1 = pcolor(struct_2D);
set(h1,'EdgeColor','none')
title(['Structured ' plot_name ])
%xlim([0 1000])
subplot(1,3,2)
h2 = pcolor(random_2D);
set(h2,'EdgeColor','none')
title(['Random ' plot_name ])
%xlim([0 1000])

%nbins = 20; %origina;
%nbins = 50; %keiko
%[c_struct, v_struct] = hist(double(struct.times), nbins);
%[c_rand, v_random] = hist(double(random.times), nbins);
bin_step = 50;
nbins = 0:bin_step:(hard_max_time-hard_min_time);
[c_struct, v_struct] = hist(double(struct.times-hard_min_time), nbins);
[c_rand, v_random] = hist(double(random.times-hard_min_time), nbins);

subplot(1,3,3)
%Q = double(struct.times-mint); %original 
%R = double(random.times-mint); %original
Q = double(struct.times-hard_min_time); %keiko
R = double(random.times-hard_min_time); %keiko

%histogram(Q) %original
histogram(Q,nbins) %keiko
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
%histogram(R) %original
histogram(R,nbins) %keiko
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
%xlim([0 2000]) %original

legend([h h1], {sprintf('Structured {%2.2f (%2.2f)}',mean(c_struct), std(c_struct)), sprintf('Random {%2.2f (%2.2f)}', mean(c_rand), std(c_rand))})



% hard_min_time = 500;
% hard_max_time = 2500;

hard_min_time = 0;
hard_max_time = 0;

% data_dir = '/home/leonardo/projects/nsdm/hill_tononi_synthesis/data';

% random_data = '/sim_1_discardtime_100.00_simtime_100.00_lambda_dg_-1.00_f_dg_20.00_retDC_30.00_retAC_30.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';
% structured_data = '/sim_1_discardtime_100.00_simtime_100.00_lambda_dg_2.00_f_dg_20.00_retDC_30.00_retAC_30.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';

% random_data = '/sim_1_discardtime_0.00_simtime_5000.00_lambda_dg_-1.00_f_dg_20.00_retDC_30.00_retAC_30.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';
% structured_data = '/sim_1_discardtime_0.00_simtime_5000.00_lambda_dg_2.00_f_dg_20.00_retDC_30.00_retAC_30.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';

% random_data = '/sim_1_discardtime_0.00_simtime_5000.00_lambda_dg_-1.00_f_dg_0.50_retDC_10.00_retAC_10.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';
% structured_data = '/sim_1_discardtime_0.00_simtime_5000.00_lambda_dg_2.00_f_dg_0.50_retDC_10.00_retAC_10.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';

data_dir = '/home/leonardo/projects/nsdm/hill_tononi_synthesis/data/erik_gvals';

% random_data = '/sim_1_discardtime_0.00_simtime_5000.00_lambda_dg_-1.00_f_dg_0.50_retDC_10.00_retAC_10.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';
% structured_data = '/sim_1_discardtime_0.00_simtime_5000.00_lambda_dg_2.00_f_dg_0.50_retDC_10.00_retAC_10.00_sim_interval_1.00_visSize_8.00_N_40.00_bar_0.00_phi_dg_0.00_detectors';

random_data = '/long_nmda_random';
structured_data = '/long_nmda_structured';

% random_data = '/ht_model_random';
% structured_data = '/ht_model_structured';

% population_name = 'Retina'; plot_name = population_name;
populations = {'Vp_v L4pyr', 'Vp_h L4pyr'}; plot_name = 'L4pyr';

struct = load([data_dir structured_data '/spikes_' populations{1} '.mat']);
struct.senders = struct.senders - min(struct.senders) + 1;
random = load([data_dir random_data '/spikes_' populations{1} '.mat']);
random.senders = random.senders - min(random.senders) + 1;

for population_name = populations(2:end)

    file_to_load = ['/spikes_' population_name{1}];

    s = load([data_dir structured_data file_to_load '.mat']);
    s.senders = s.senders - min(s.senders) + 1;

    struct.senders = [struct.senders  s.senders+max(struct.senders)+1];
    struct.times = [struct.times s.times];

    r = load([data_dir random_data file_to_load '.mat']);
    r.senders = r.senders - min(r.senders) + 1;
    
    random.senders = [random.senders r.senders+max(random.senders)+1];
    random.times = [random.times r.times];
    
end

time_grain = 1; % mili seconds
% time_grain = 1000; %micro senconds

%convert time scale
struct.times = int32(time_grain * struct.times);
random.times = int32(time_grain * random.times);

% apply threshold
if hard_max_time ~= 0
    too_late = struct.times < hard_max_time;
    struct.times(too_late) = [];
    struct.senders(too_late) = [];
    
    too_late = random.times < hard_max_time;
    random.times(too_late) = [];
    random.senders(too_late) = [];
end

if hard_min_time ~= 0
    too_soon = struct.times > hard_min_time;
    struct.times(too_soon) = [];
    struct.senders(too_soon) = [];
    
    too_soon = random.times > hard_min_time;
    random.times(too_soon) = [];
    random.senders(too_soon) = [];
end


maxt = int32((100*time_grain)*(ceil(max(struct.times)/(100*time_grain))));
if maxt ~= int32((100*time_grain)*(ceil(max(random.times)/(100*time_grain))))
    error('random and structured do not have the same duration!')
end

% discard first part of simulation
mint = int32(time_grain*min(struct.times));
if mint ~= int32(time_grain*min(random.times))
    error('random and structured do not have the same startnig point!')
end

maxn = int32(100*(ceil(max(struct.senders)/100)));
if maxn ~= int32(100*(ceil(max(random.senders)/100)))
    error('random and structured do not have the same # of neurons!')
end

minn = int32(min(struct.senders));
if minn ~= int32(min(random.senders))
    error('random and structured do not have the same startnig point!')
end

struct.senders = int32(struct.senders);
random.senders = int32(random.senders);

fprintf('\n\nEstimated memory needed : %2.2f Gb\n\n', (2*maxt*maxn*8) / (2^30))


%% plot


% downsample to visualitaion (screen resolution)

sizeN = 1600;
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


figure


subplot(1,3,1)
h1 = pcolor(struct_2D);
set(h1,'EdgeColor','none')
title(['Structured ' plot_name ])
xlim([0 1000])
subplot(1,3,2)
h2 = pcolor(random_2D);
set(h2,'EdgeColor','none')
title(['Random ' plot_name ])
xlim([0 1000])

nbins = 20;
[c_struct, v_struct] = hist(double(struct.times), nbins);
[c_rand, v_random] = hist(double(random.times), nbins);

subplot(1,3,3)
Q = double(struct.times-mint);
R = double(random.times-mint);
histogram(Q)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
histogram(R)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
xlim([0 2000])

legend([h h1], {sprintf('Structured {%2.2f (%2.2f)}',mean(c_struct), std(c_struct)), sprintf('Random {%2.2f (%2.2f)}', mean(c_rand), std(c_rand))})


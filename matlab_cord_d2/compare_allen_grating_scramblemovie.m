clear all

dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/allen_tmp/';

grating_filename = 'grat_7_0_spks.dat';
%grating_filename = 'mov_2_0_spks.dat';
%scramble_filename = 'mov_2_scrbl_xy_0_spks.dat';
scramble_filename = 'mov_2_scrbl_t_0_spks.dat';

% Set time window for analysis based on Yazan's email
%hard_min_time = 500.0;
%hard_max_time = 3000.0;
hard_min_time = 0.0;
hard_max_time = 3000.0;

% Set neuron index for analysis based on Yazan's email
% Here I use only realistic excitatory neurons
populations_idx = 0:1:45000; 
plot_name = 'All neurons';
%
%populations_idx = 0:1:8499; 
%plot_name = 'Realistic Exc';
%
%populations_idx = 8500:1:9999; 
%plot_name = 'PV-positive (Inh)';
%
%populations_idx = 10000:1:39749; 
%plot_name = 'LIF (Exc)';
%
%populations_idx = 39750:1:45000; 
%plot_name = 'LIF (Inh)';

% Load data
tmp_grating = load( strcat(dirname,grating_filename) );
tmp_scramble = load( strcat(dirname,scramble_filename) );

% Make structure
% This is not necessory, but I prefer this to increase readability of this code.
grating = struct('times', tmp_grating(:,1), 'senders', tmp_grating(:,2));
scramble = struct('times', tmp_scramble(:,1), 'senders', tmp_scramble(:,2));

% Convert neuon index
% Allen format starts from 0, matlab starts from 1
grating.senders = grating.senders + 1;
scramble.senders = scramble.senders + 1;

%clear tmp_grating
%clear tmp_scramble

% Select data within time window (hard_min_time < t < hard_max_time)
if hard_max_time ~= 0
    too_late = grating.times > hard_max_time;
    grating.times(too_late) = [];
    grating.senders(too_late) = [];
    
    too_late = scramble.times > hard_max_time;
    scramble.times(too_late) = [];
    scramble.senders(too_late) = [];
end

if hard_min_time ~= 0
    too_soon = grating.times < hard_min_time;
    grating.times(too_soon) = [];
    grating.senders(too_soon) = [];
    
    too_soon = scramble.times < hard_min_time;
    scramble.times(too_soon) = [];
    scramble.senders(too_soon) = [];
end

% Select target neurons
tgt_rows_grating = ismember(grating.senders, populations_idx);
grating.times = grating.times(tgt_rows_grating);
grating.senders = grating.senders(tgt_rows_grating);

tgt_rows_scramble = ismember(scramble.senders, populations_idx);
scramble.times = scramble.times(tgt_rows_scramble);
scramble.senders = scramble.senders(tgt_rows_scramble);

%%
time_grain = 1; %mili seconds
%time_grain = 1000; %micro senconds

%convert time scale
grating.times = int32(time_grain * grating.times);
scramble.times = int32(time_grain * scramble.times);

%{
maxt = int32((100*time_grain)*(ceil(max(grating.times)/(100*time_grain))));
if maxt ~= int32((100*time_grain)*(ceil(max(scramble.times)/(100*time_grain))))
    error('random and structured do not have the same duration!')
end

% discard first part of simulation
mint = int32(time_grain*min(grating.times));
if mint ~= int32(time_grain*min(scramble.times))
    error('random and structured do not have the same startnig point!')
end

maxn = int32(100*(ceil(max(grating.senders)/100)));
if maxn ~= int32(100*(ceil(max(scramble.senders)/100)))
    error('random and structured do not have the same # of neurons!')
end

minn = int32(min(grating.senders));
if minn ~= int32(min(scramble.senders))
    error('random and structured do not have the same startnig point!')
end
%}
mint = hard_min_time;
maxt = hard_max_time;
minn = min(populations_idx);
maxn = max(populations_idx);

grating.senders = int32(grating.senders);
scramble.senders = int32(scramble.senders);

fprintf('\n\nEstimated memory needed : %2.2f Gb\n\n', (2*maxt*maxn*8) / (2^30))


%% plot
N_neuron = length(populations_idx);
N_sample = hard_max_time - hard_min_time + 1;
grating_2D_all = zeros(N_neuron, N_sample);
scramble_2D_all = zeros(N_neuron, N_sample);

column = 1;
for t = hard_min_time:1:hard_max_time
    tmp_grating_rows = find(grating.times==t);
    tmp_grating_senders = grating.senders(tmp_grating_rows);
    grating_2D_all(tmp_grating_senders,column) = 1;

    tmp_scramble_rows = find(scramble.times==t);
    tmp_scramble_senders = scramble.senders(tmp_scramble_rows);
    scramble_2D_all(tmp_scramble_senders,column) = 1;

    column = column + 1;
end

%draw_neuron_idx = 1:15:N_neuron;
%draw_time_idx = 1:2:N_sample;
draw_neuron_idx = 1:50:N_neuron;
draw_time_idx = 1:2:N_sample;
grating2D = grating_2D_all(draw_neuron_idx,draw_time_idx);
scramble2D = scramble_2D_all(draw_neuron_idx,draw_time_idx);

%{
% downsample to visualitaion (screen resolution)

%sizeN = 1600; %original
%sizeT = 2001; %original
sizeN = 1000;
sizeT = 2001;

downsampleN = ceil(double(maxn-minn+1)/double(sizeN)); 
downsampleT = ceil(double(maxt-mint+1)/double(sizeT)); 

struct_2D = zeros(sizeN, sizeT);
random_2D = struct_2D;

for is = 1:length(grating.senders)
    struct_2D(ceil(double(grating.senders(is)-minn+1)/double(downsampleN)), ceil(double(grating.times(is)-mint+1)/double(downsampleT))) = 1; 
end

for is = 1:length(scramble.senders)
    random_2D(ceil(double(scramble.senders(is)-minn+1)/double(downsampleN)), ceil(double(scramble.times(is)-mint+1)/double(downsampleT))) = 1; 
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

%{
subplot(1,3,1)
h1 = pcolor(grating2D);
set(h1,'EdgeColor','none')
title(['Structured ' plot_name ])
%xlim([0 1000])
subplot(1,3,2)
h2 = pcolor(scramble2D);
set(h2,'EdgeColor','none')
title(['Random ' plot_name ])
%xlim([0 1000])
%}

%nbins = 20; %origina;
bin_step = 50;
nbins = 0:bin_step:(hard_max_time-hard_min_time);
[c_grating, v_grating] = hist(double(grating.times-hard_min_time), nbins);
[c_scramble, v_scramble] = hist(double(scramble.times-hard_min_time), nbins);


%subplot(1,3,3)
%Q = double(grating.times-mint);
%R = double(scramble.times-mint);
Q = double(grating.times-hard_min_time); %keiko
R = double(scramble.times-hard_min_time); %keiko

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

legend([h h1], {sprintf('Structured {%2.2f (%2.2f)}',mean(c_grating), std(c_grating)), sprintf('Random {%2.2f (%2.2f)}', mean(c_scramble), std(c_scramble))})
title(plot_name)


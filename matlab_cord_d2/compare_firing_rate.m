clear

hard_min_time = 500.0;
%hard_min_time = 0.0;
hard_max_time = 3500.0;

%hard_min_time = 0;
%hard_max_time = 0;

% data_dir = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency';
% random_data = '/random_rate150';
% structured_data = '/2Hz_vertical_rate150';
%
data_dir = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern';
random_data = '/walking_human_rate100_scrbl_xy';
%structured_data = '/moving_rectangle2_rate100_dst';
%random_data = '/walking_human_rate20_scrbl_xy';
structured_data = '/walking_human_rate100_dst';

% --- Load vertical data (struct/random)
struct_v = load([data_dir, structured_data, '/spike_', 'Vp_v_L4_exc', '.mat']);
struct_v.senders = struct_v.senders - min(struct_v.senders) + 1;
random_v = load([data_dir, random_data, '/spike_', 'Vp_v_L4_exc', '.mat']);
random_v.senders = random_v.senders - min(random_v.senders) + 1;

% --- Load horizontal data (struct/random)
struct_h = load([data_dir, structured_data, '/spike_', 'Vp_h_L4_exc', '.mat']);
struct_h.senders = struct_h.senders - min(struct_h.senders) + 1;
random_h = load([data_dir, random_data, '/spike_', 'Vp_h_L4_exc', '.mat']);
random_h.senders = random_h.senders - min(random_h.senders) + 1;

% --- Set all data
struct_all.senders = [struct_v.senders struct_h.senders];
random_all.senders = [random_v.senders random_h.senders];
struct_all.times = [struct_v.times struct_h.times];
random_all.times = [random_v.times random_h.times];


time_grain = 1; % mili seconds
% time_grain = 1000; %micro senconds

%convert time scale
struct_v.times = int32(time_grain * struct_v.times);
struct_h.times = int32(time_grain * struct_h.times);
struct_all.times = int32(time_grain * struct_all.times);

random_v.times = int32(time_grain * random_v.times);
random_h.times = int32(time_grain * random_h.times);
random_all.times = int32(time_grain * random_all.times);

% apply threshold
if hard_max_time ~= 0
    %--- Structured (v/h/all)
    too_late = struct_v.times > hard_max_time;
    struct_v.times(too_late) = [];
    struct_v.senders(too_late) = [];
    
    too_late = struct_h.times > hard_max_time;    
    struct_h.times(too_late) = [];
    struct_h.senders(too_late) = [];
    
    too_late = struct_all.times > hard_max_time;        
    struct_all.times(too_late) = [];
    struct_all.senders(too_late) = [];
    
    %--- random (v/h/all)
    too_late = random_v.times > hard_max_time;
    random_v.times(too_late) = [];
    random_v.senders(too_late) = [];

    too_late = random_h.times > hard_max_time;
    random_h.times(too_late) = [];
    random_h.senders(too_late) = [];

    too_late = random_all.times > hard_max_time;
    random_all.times(too_late) = [];
    random_all.senders(too_late) = [];
end

if hard_min_time ~= 0
    
    %--- Structured (v/h/all)
    too_soon = struct_v.times < hard_min_time;
    struct_v.times(too_soon) = [];
    struct_v.senders(too_soon) = [];
    
    too_soon = struct_h.times < hard_min_time;
    struct_h.times(too_soon) = [];
    struct_h.senders(too_soon) = [];
    
    too_soon = struct_all.times < hard_min_time;
    struct_all.times(too_soon) = [];
    struct_all.senders(too_soon) = [];
    
    %--- Random (v/h/all)
    too_soon = random_v.times < hard_min_time;
    random_v.times(too_soon) = [];
    random_v.senders(too_soon) = [];

    too_soon = random_h.times < hard_min_time;
    random_h.times(too_soon) = [];
    random_h.senders(too_soon) = [];

    too_soon = random_all.times < hard_min_time;
    random_all.times(too_soon) = [];
    random_all.senders(too_soon) = [];

end


%% plot
figure

bin_step = 50;
nbins = 0:bin_step:(hard_max_time-hard_min_time);

%--- All population (vertical + horizontal)
subplot(1,3,1)
[c_struct_all, v_struct_all] = hist(double(struct_all.times-hard_min_time), nbins);
[c_random_all, v_random_all] = hist(double(random_all.times-hard_min_time), nbins);

Q_all = double(struct_all.times-hard_min_time);
R_all = double(random_all.times-hard_min_time);
tmpQ = histogram(Q_all,nbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
tmpR = histogram(R_all,nbins);
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
title('All popurations')
legend([h h1], {sprintf('Structured {%2.2f (%2.2f)}',...
                mean(c_struct_all), std(c_struct_all)),...
                sprintf('Random {%2.2f (%2.2f)}',...
                mean(c_random_all), std(c_random_all))})

max_val = max( max(tmpQ.Values), max(tmpR.Values) );
max_axis = 10 * ceil(0.1*max_val);
ylim([0 max_axis])
            
%--- vertical only
subplot(1,3,2)
[c_struct_v, v_struct_v] = hist(double(struct_v.times-hard_min_time), nbins);
[c_random_v, v_random_v] = hist(double(random_v.times-hard_min_time), nbins);

Q_v = double(struct_v.times-hard_min_time);
R_v = double(random_v.times-hard_min_time);
histogram(Q_v,nbins)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
histogram(R_v,nbins)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
title('Vertical popuration')
legend([h h1], {sprintf('Structured {%2.2f (%2.2f)}',...
                mean(c_struct_v), std(c_struct_v)),...
                sprintf('Random {%2.2f (%2.2f)}',...
                mean(c_random_v), std(c_random_v))})
ylim([0 max_axis])

%--- horizontal only
subplot(1,3,3)
[c_struct_h, v_struct_h] = hist(double(struct_h.times-hard_min_time), nbins);
[c_random_h, v_random_h] = hist(double(random_h.times-hard_min_time), nbins);

Q_h = double(struct_h.times-hard_min_time);
R_h = double(random_h.times-hard_min_time);
histogram(Q_h,nbins)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
histogram(R_h,nbins)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.75);
title('Horizontal popuration')
legend([h h1], {sprintf('Structured {%2.2f (%2.2f)}',...
                mean(c_struct_v), std(c_struct_v)),...
                sprintf('Random {%2.2f (%2.2f)}',...
                mean(c_random_h), std(c_random_h))})
ylim([0 max_axis])
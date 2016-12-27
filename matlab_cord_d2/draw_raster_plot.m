clear

%exp_dir = '2Hz_vertical_rate100';
exp_dir = 'random_rate100';
%exp_dir = 'random_rate100_alldata';
%exp_dir = 'moving_rectangle2_rate100_dst';
%exp_dir = 'walking_human_rate100_dst';
%exp_dir = 'walking_human_rate100_scrbl_xy';

%
%filename = strcat('/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/',...
filename = strcat('/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/',...
                   exp_dir,...
                  '/spike_Vp_v_L4_exc.mat');
src_data = load(filename);

time_threshold_min = 500;
time_threshold_max = 3500;

% Remove unnecessory part of the data
early_idx = find(src_data.times<time_threshold_min);
src_data.senders(early_idx) = [];
src_data.times(early_idx) = [];

late_idx = find(src_data.times>time_threshold_max);
src_data.senders(late_idx) = [];
src_data.times(late_idx) = [];

% Convert index from NEST to MATLAB
src_data.senders = src_data.senders - min(src_data.senders) + 1;
src_data.times = src_data.times - min(src_data.times) + 1;

% Put spike data into a matrix
num_neuron = max(src_data.senders) - min(src_data.senders) + 1;
duration = max(src_data.times) - min(src_data.times);
spike_mat = zeros(num_neuron,duration);
for t = 1:1:duration
    tmp_t_idx = find(src_data.times==t);
    tmp_n_idx = src_data.senders(tmp_t_idx);
    spike_mat(tmp_n_idx,t) = 1;
end

%subplot(1,2,1)
%plot(sum(spike_mat,1))

%subplot(1,2,2)
figure
scatter(src_data.times, src_data.senders, '.')

%--- Frequency analysis
X = sum(spike_mat,1);

Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 3000;             % Length of signal
t = (0:L-1)*T;        % Time vector
%plot(1000*t(1:50),X(1:50))
%title('Signal Corrupted with Zero-Mean Random Noise')
%xlabel('t (milliseconds)')
%ylabel('X(t)')

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure
%plot(f,P1)
plot(f(2:end),P1(2:end))
title(strrep(exp_dir,'_',' '))
xlabel('f (Hz)')
ylabel('|P1(f)|')
grid on


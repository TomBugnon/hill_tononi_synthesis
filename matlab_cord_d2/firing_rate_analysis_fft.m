clear

%% Set parameters
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/2Hz_vertical_rate100/';
dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_frequency/random_rate100/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human_rate100_dst/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human_rate100_scrbl_xy/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/walking_human_rate100_scrbl_t/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/moving_rectangle2_rate100_dst/';
%dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/lobustness_pattern/scramble_test_moving_rectangle2_rate20_dst/';

filename_v = strcat(dirname, 'spike_Vp_v_L4_exc.mat');
filename_h = strcat(dirname, 'spike_Vp_h_L4_exc.mat');

time_threshold_min = 500;
time_threshold_max = 5000;

time_window_min = 1; %[ms]
time_window_max = 300; %[ms]
window_step = 2;


%% Load data
% Output from NEST uses sparse-expression, so I converted to matrix
% spike_mat
% --- row: neuron index
% --- col: time sequence
% --- If i-th neuron fired at time j, spike_mat(i,j)==1, otherwise 0
% --- Note that time is limited from "time_window_min" to "time_window_max"
% --- so, j-th column in spike_mat represents the time of (j + time_window_min) in NEST simulator.
spike_mat_v = convert_nestspike2mat(filename_v,time_threshold_min,time_threshold_max);
spike_mat_h = convert_nestspike2mat(filename_h,time_threshold_min,time_threshold_max);

spike_mat = [spike_mat_v; spike_mat_h];

%% FFT
X = sum(spike_mat,1);

Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = size(X,2);        % Length of signal
t = (0:L-1)*T;        % Time vector
%plot(1000*t(1:50),X(1:50))
%title('Signal Corrupted with Zero-Mean Random Noise')
%xlabel('t (milliseconds)')
%ylabel('X(t)')

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
f = Fs*(0:(L/2))/L;

figure
plot(f(2:end),P1(2:end))
xlabel('f (Hz)')
ylabel('|P1(f)|')
grid on





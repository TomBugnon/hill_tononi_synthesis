

clear

%% Set parameters

dirname = '/Users/Tom/Desktop/test_leonardo_network/network_full_leonardo_intact_edge_wrap_1_Np_40_Ns_30_p_ratio_2/vertical_rate100_run1/'
mintime = 500 ;
maxtime = 2500 ;

%Which ones are we interested in, and % How do we plot them (color for layer, dotted/full line for Vs/Vp
pops = {'Vp_v_L23_exc', 'Vp_v_L4_exc', 'Vp_v_L56_exc', 'Vs_v_L23_exc', 'Vs_v_L4_exc', 'Vs_v_L56_exc'}
linestyle = {'-', '-', '-', '--', '--', '--'}
pops = {'Vs_v_L23_exc', 'Vs_v_L4_exc', 'Vs_v_L56_exc'}
linestyle = {'-', '-', '-', '--', '--', '--'}


%% load
pops_filenames = arrayfun(@(popname)(strcat(dirname, 'spike_', popname, '.mat')), pops)

pops_activity = cellfun(@(filename)(convert_nestspike2mat(filename,mintime,maxtime)), pops_filenames, 'UniformOutput', false)




%% Plot 

fig = figure;
hold on
sub1 = subplot( 2, 1, 1)
sub2 = subplot( 2, 1, 2)


%% FFT & TC


Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period


hfft = []
htc = []
for popnumb = 1:numel(pops_activity)
    spike_mat = pops_activity{popnumb};
    

    X = sum(spike_mat,1);
    
    Xsmooth = smooth(X)
    
  
    L = size(X,2);        % Length of signal
    t = (0:L-1)*T;        % Time vector

    Y = fft(X);
    P2 = abs(Y/L);
    P1 = smooth(P2(1:L/2+1));
    f = Fs*(0:(L/2))/L;
    
    axes(sub1)
    hold on
    times = ((maxtime - mintime)/L) * mintime:(maxtime-1);  %Perhaps wrong 
    htc(popnumb) = plot(times, Xsmooth, 'linewidth', 1.5, 'LineStyle', linestyle{popnumb})
    hold off
    axes(sub2)
    hold on
    hfft(popnumb) = plot(f(2:end),P1(2:end), 'linewidth', 1.5,'LineStyle', linestyle{popnumb})
    hold off
    
    
end

axes(sub1)
title('Smoothed time courses','FontSize', 25)
xlabel('time', 'FontSize', 25)
ylabel('# spikes', 'FontSize', 25)
legend(htc, pops, 'FontSize', 25)


axes(sub2)
title('Smoothed FFT decomposition', 'FontSize', 25)
xlabel('f (Hz)', 'FontSize', 25)
ylabel('|P1(f)|', 'FontSize', 25)
xlim([0,100])

legend(hfft, pops, 'FontSize', 25)


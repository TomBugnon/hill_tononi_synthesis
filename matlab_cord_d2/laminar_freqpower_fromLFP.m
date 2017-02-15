

clear

%% Set parameters

dirname = '/Users/Tom/Desktop/test_leonardo_network/network_full_leonardo_intact_edge_wrap_1_Np_40_Ns_30_p_ratio_2/vertical_rate100_run1/'
mintime = 500 ;
maxtime = 2500 ;

LFPmethod = 'Isyn' % in {'Isyn', 'Vm'}
LFPmethod = 'Vm'


%Which ones are we interested in, and % How do we plot them (color for layer, dotted/full line for Vs/Vp
pops = {'Vp_v_L23_exc', 'Vp_v_L4_exc', 'Vp_v_L56_exc', 'Vs_v_L23_exc', 'Vs_v_L4_exc', 'Vs_v_L56_exc'}
linestyle = {'-', '-', '-', '--', '--', '--'}
colors = {'red', 'blue', 'green', 'red', 'blue', 'green'}


%Which fields are we interested in loading
fields = {'V_m', 'g_GABA_A', 'g_AMPA', 'g_GABA_A', 'g_GABA_B', 'g_NMDA'}
%TODO: should import from model
reversal_potentials = struct('g_GABA_A', -70, 'g_GABA_B', -90, 'g_AMPA', 0, 'g_NMDA', 0)




%% load
pops_filenames = arrayfun(@(popname)(strcat(dirname, 'recorder_', popname, '.mat')), pops)



%convert. structure for each population, with fields for each recorded
%variable
pops_activity = cellfun(@(filename)(convert_nestrecorder2mat(filename,mintime,maxtime, fields)), pops_filenames, 'UniformOutput', false)



%% Compute LFP
pops_lfp = cell(size(pops));
for popnum = 1:numel(pops)
    lfp = zeros(size(pops_activity{popnum}.V_m));
    synapses = fieldnames(reversal_potentials);
    
    
    if strcmp(LFPmethod, 'Isyn')
        
        for synapsenumb = 1:numel(synapses)
            syn = synapses{synapsenumb};
            %add sum of absolute synaptic currents (conductance * (vm -
            %reversal potential)
            conductance = pops_activity{popnum}.(syn);
            erev = reversal_potentials.(syn);
            vm = pops_activity{popnum}.V_m;
            lfp = lfp + abs(conductance .* (vm - erev));
            
            
        end
    else 
        lfp = pops_activity{popnum}.V_m;

    end
        
        
    pops_lfp{popnum} = lfp;
    
end    


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
for popnumb = 1:numel(pops_lfp)
    lfp_matrix = pops_lfp{popnumb};
    

    X = mean(lfp_matrix,1);
    
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
    htc(popnumb) = plot(times, X, 'linewidth', 1.5, 'LineStyle', linestyle{popnumb}, 'color', colors{popnumb})
    xlim([mintime, maxtime])
    hold off
    axes(sub2)
    hold on
    hfft(popnumb) = plot(f(2:end),P1(2:end), 'linewidth', 1.5,'LineStyle', linestyle{popnumb}, 'color', colors{popnumb})
    hold off
    
    
end

axes(sub1)
title('Smoothed time courses','FontSize', 25)
xlabel('time', 'FontSize', 25)
ylabel('LFP (Sum of absolute synaptic currents)', 'FontSize', 25)
legend(htc, pops, 'FontSize', 25)


axes(sub2)
title('Smoothed FFT decomposition', 'FontSize', 25)
xlabel('f (Hz)', 'FontSize', 25)
ylabel('|P1(f)|', 'FontSize', 25)
xlim([0,100])

legend(hfft, pops, 'FontSize', 25)



dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/keep/';
%dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/examine_EPSP_IPSP_keiko_nov03_pm03/';
%dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/examine_EPSP_IPSP_pablo/';
filename = 'recorder_Vp_v_L23_exc.mat';
%filename = 'recorder_Vp_v_L4_exc.mat';
%filename = 'recorder_Vp_v_L56_exc.mat';
data = load( strcat(dirname,filename) );

min_neuron_idx = min(data.senders);
max_neuron_idx = max(data.senders);

num_neurons = max_neuron_idx - min_neuron_idx + 1;
t_length = length(find(data.senders==min_neuron_idx));
pos = zeros(num_neurons, t_length);
I_syn_NMDA = zeros(num_neurons, t_length);
I_syn_AMPA = zeros(num_neurons, t_length);
I_syn_GABAA = zeros(num_neurons, t_length);
I_syn_GABAB = zeros(num_neurons, t_length);
g_AMPA = zeros(num_neurons, t_length);
g_NMDA = zeros(num_neurons, t_length);
g_GABAA = zeros(num_neurons, t_length);
g_GABAB = zeros(num_neurons, t_length);
V_m = zeros(num_neurons, t_length);

for i = 1:1:num_neurons
    
    target_idx = min_neuron_idx + i -1;
    tmp_pos = find(data.senders==target_idx);
    pos(i,:) = tmp_pos;
    
    I_syn_NMDA(i,:) = data.I_syn_NMDA(tmp_pos);
    I_syn_AMPA(i,:) = data.I_syn_AMPA(tmp_pos);
    I_syn_GABAA(i,:) = data.I_syn_GABA_A(tmp_pos);
    I_syn_GABAB(i,:) = data.I_syn_GABA_B(tmp_pos);    
    g_AMPA(i,:)  = data.g_AMPA(tmp_pos);
    g_NMDA(i,:)  = data.g_NMDA(tmp_pos);
    g_GABAA(i,:) = data.g_GABAA(tmp_pos);
    g_GABAB(i,:) = data.g_GABAB(tmp_pos);
    
    V_m(i,:) = data.V_m(tmp_pos);

end


I_syn_thredhold = 1.0;
g_thredhold = 0.05;
na_ratio = zeros(num_neurons, t_length);
ei_ratio = zeros(num_neurons, t_length);
for i = 1:1:num_neurons
    na_tmp_idx = find(I_syn_AMPA(i,:)>I_syn_thredhold);
    na_ratio(i,na_tmp_idx) = I_syn_NMDA(i,na_tmp_idx) ./ I_syn_AMPA(i,na_tmp_idx);
    
    ei_tmp_idx = find(g_AMPA(i,:)>g_thredhold);
    ei_ratio(i,ei_tmp_idx) = g_GABAA(i,ei_tmp_idx) ./ g_AMPA(i,ei_tmp_idx);
    %exc_data = g_AMPA + g_NMDA;
    %inh_data = g_GABAA + g_GABAB;
    %ei_tmp_idx = find(exc_data(i,:)>g_thredhold);
    %ei_ratio(i,ei_tmp_idx) = inh_data(i,ei_tmp_idx) ./ exc_data(i,ei_tmp_idx);

end

% draw figure (NMDA-AMPA ratio)
na_mean = mean(na_ratio);
figure; plot(na_mean)
grid on
title('NMDA-AMPA ratio in pablo full network')
xlabel('Time[ms]')
ylabel('I_N_M_D_A / I_A_M_P_A')

% draw figure (EI ratio)
ei_mean = mean(ei_ratio);
figure; plot(ei_mean)
grid on
title('Excitatory-Inhibitory ratio in pablo full network')
xlabel('Time[ms]')
ylabel('g_G_A_B_A_A / g_A_M_P_A')


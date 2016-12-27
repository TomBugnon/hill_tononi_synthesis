%dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/examine_EPSP_IPSP_keiko/';
%dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/examine_EPSP_IPSP_keiko_nov03_pm03/';
dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/';
%target_gid = 14406+2; % Vp_v_L4_exc  %pablo
%target_area = 'Vp_v_L4_exc';

target_gid = 25606+20; % Vp_v_L4_exc %keiko
target_area = 'Vp_v_L4_exc';

%target_gid = 30406+20; % Vp_v_L56_exc %keiko
%target_area = 'Vp_v_L56_exc';

%target_gid = 20806+20; % Vp_v_L56_exc %keiko
%target_area = 'Vp_v_L23_exc';

%synapse = 'AMPA';
%synapse = 'GABA_A';
synapse = 'GABA_B';
population = {'Vp_v', 'Vp_h', 'Tp','Vs_v', 'Vs_h'};
model = {'L23_exc', 'L23_inh',...
         'L4_exc', 'L4_inh',...
         'L56_exc', 'L56_inh',...
         'Tp_exc', 'Tp_inh'};

% Find globalIDs of pre_neurons
pre_all = [];
for p = 1:1:length(population)
    connection_fname = strcat(dirname, 'connection_', population{p}, '_', synapse, '_syn.dat');
    if exist(connection_fname,'file')>0
        pop = load(connection_fname);
    end
    if ~isempty(pop)
        rows = find(pop(:,2)==target_gid);
        pre = pop(rows,1);
        pre_all = [pre_all; pre];
    end
end

% Set conductance values
conductance_all = [];
for p = 1:1:length(population)
    
    for m = 1:1:length(model)
        
        % Load data    
        conductance_fname = strcat(dirname, 'recorder_', population{p}, '_', model{m}, '.mat');   
        
        if exist(conductance_fname,'file')>0

            data = load(conductance_fname);
            
            % Preparation
            min_neuron_idx = min(data.senders);
            max_neuron_idx = max(data.senders);
            num_neurons = max_neuron_idx - min_neuron_idx + 1;
            t_length = length(find(data.senders==min_neuron_idx));        
            conductance = zeros(num_neurons, t_length);   
                
            % Get conductance data
            % 'conductance' matrix's row represents neuron index(local-ID)
            % 'conductance' matrix's col represents time seaquence    
            for i = 1:1:num_neurons

                target_idx = min_neuron_idx + i -1;
                tmp_pos = find(data.senders==target_idx);

                if strcmp(synapse,'AMPA')
                    conductance(i,:)  = data.g_AMPA(tmp_pos);
                elseif strcmp(synapse,'NMDA')
                    conductance(i,:)  = data.g_NMDA(tmp_pos);            
                elseif strcmp(synapse,'GABA_A')
                    conductance(i,:)  = data.g_GABAA(tmp_pos);            
                elseif strcmp(synapse,'GABA_B')
                    conductance(i,:)  = data.g_GABAB(tmp_pos);
                else
                    disp('Specify synapse type')
                end

                %if strcmp( target_area, strcat(population{p},'_', model{m}))
                if target_idx == target_gid
                    V_m = data.V_m(tmp_pos);
                end
            end
            
            %{
            % If the target(post) neuron is in the area processing now, get V_m data
            if strcmp( target_area, strcat(population{p},'_', model{m}))
                strcat(population{p},'_', model{m})
                target_pos = find(data.senders==target_idx);
                V_m = data.V_m(target_pos);
            end
            %}
            
            % Find pre_neuron position in 'conductance' matrix
            tmp = ismember(data.senders, pre_all);
            if sum(tmp)>0
                
                conductance_fname                
                pre_in_roi_gid = unique( data.senders(tmp) );
                pre_in_roi_lid = pre_in_roi_gid - min_neuron_idx + 1;
                length(pre_in_roi_lid)
                
                %{
                % Draw time seaquence of pre-neuron's conductances
                figure
                for i=1:1:length(pre_in_roi_lid)
                    plot( conductance(pre_in_roi_lid(i),:) );
                    hold on
                end
                title_str = strrep( strcat(population{p}, '-', model{m}, '-', synapse), '_', '-');
                title(title_str)
                grid on
                hold off
                %}
                total_conductance_area = sum(conductance,1);
                conductance_all = [conductance_all; total_conductance_area];
            
            end            
        end
    end
end

figure
plot(V_m)
grid on
title( 'V_m of the target(post) neuron')

% Calculate post synaptic potential (EPSP/IPSP)
if strcmp(synapse,'AMPA')
    Erev = 0;
elseif strcmp(synapse,'NMDA')
    Erev = 0;
elseif strcmp(synapse,'GABA_A')
    Erev = -70;
elseif strcmp(synapse,'GABA_B')
    Erev = -90;
else
    disp('Specify synapse type')
end

PSP = -1 * conductance_all .* repmat((V_m - Erev),size(conductance_all,1),1);

figure; 
plot(PSP')
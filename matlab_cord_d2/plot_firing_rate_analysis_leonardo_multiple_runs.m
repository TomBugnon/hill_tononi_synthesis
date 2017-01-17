
load([root_dir '/last_group_analyses.mat'], 'all_x_time', 'all_max_resolution', 'all_max_d2',  'all_mean_p', 'all_std_p',... 
                                            'all_threshold_all', 'all_max_threshold');

%% Display results

fs = cell(1,ntemplates);
for id = 1:ntemplates
    fs{id} = figure;
    set(gcf, 'Name', dirnames_template{id})
end

avg_max_d2 = zeros(id,run);
avg_max_resolution = zeros(id,run);

for run = 1:nruns
    
    fprintf('Plotting run %d of %d\n', run, nruns)
    for id = 1:ntemplates
        
        figure(fs{id})

        x_time = all_x_time{id, run};
        max_resolution = all_max_resolution{id, run};
        max_d2 = all_max_d2{id, run};
        mean_p = all_mean_p{id, run};
        std_p = all_std_p{id, run};
        
        avg_max_d2(id,run) = max_d2;
        avg_max_resolution(id,run) = max_resolution;

        threshold_all = all_threshold_all{id, run};
        max_threshold = all_max_threshold{id, run};

        if fixed_threshold
            sublayout = [121, 122, 0];
        else
            sublayout = [131, 132, 133];
        end

        subplot(sublayout(1))
        hold on
        plot(x_time, d2_all)
        hold on; 
        plot(max_resolution, max_d2, 'ob', 'MarkerSize',10,...
        'MarkerEdgeColor','k','MarkerFaceColor','g')
        grid on
        title( strcat('D2 (maxD2=',num2str(max_d2),')') );
        set(gca,'XTick', x_time(1:80:end))
        %set(gca,'XTickLabel', x_str(1:20:end))
        %ylim([0 6000])
        ylim([0 1])

        subplot(sublayout(2))
        hold on
        %bar(mean(max_p_positive))
        % hold on
        % errorbar(mean(max_p), std(max_p))
        % title( strcat('max probability (timewindow = ',num2str(max_resolution), ')') )
        %plot(x_val, mean_p(:,1), 'k'); 
        silent_h = shadedErrorBar(x_time, mean_p(:,1),std_p(:,1),'k',1);
        hold on; 
        %plot(x_val, mean_p(:,2), 'r')
        active_h = shadedErrorBar(x_time, mean_p(:,2),std_p(:,2),'r',1);
        legend([silent_h.mainLine, active_h.mainLine], {'silent', 'active'})
        plot(max_resolution*[1, 1], mean_p(x_time==max_resolution,:), 'ob', 'MarkerSize',10,...
        'MarkerEdgeColor','k','MarkerFaceColor','g')
    %     plot(max_resolution*[1, 1], mean(max_p_positive), 'ob', 'MarkerSize',10,...
    %     'MarkerEdgeColor','k','MarkerFaceColor','g')
        set(gca,'XTick', x_time(1:80:end))
        ylim([0 1.0])
        title( strcat('probability (win maxD2 = ',num2str(max_resolution), ')') )

        if sublayout(3)
            subplot(1,3,3)
            plot(x_time,threshold_all)
            hold on; 
            plot(max_resolution, max_threshold, 'ob', 'MarkerSize',10,...
            'MarkerEdgeColor','k','MarkerFaceColor','g')
            grid on
            title( strcat('Threshold value (firing rate = ',num2str(max_threshold),')') )
            set(gca,'XTick', x_time(1:20:end))
        %     ylim([0 25])
        end

    %     title_str = strsplit(dirname, '/');
    %     set(gcf, 'Name', title_str{end-1})

        % figure
        % histogram(max_p(:,1),20,'FaceColor','k')
        % hold on
        % histogram(max_p(:,2),20,'FaceColor','r')
        % title(max_resolution)

        %{
        figure;
        bar(firing_hist_all(10:10:end,:)')
        figure;
        bar(median_each_all(10:10:end,:)')
        %}

        %for i=10:10:size(firing_rate_all,1)
        %    bar(firing_rate_all(i,:));
        %    hold on
        %end


        % %% Draw raster plot & firing_rate
        % num_draw_t = 1000;
        % %num_draw_neuron = 3200;
        % 
        % frate = zeros(1,num_draw_t);
        % figure
        % subplot(2,1,1)
        % to_plot_x = [];
        % to_plot_y = [];
        % for t=1:1:num_draw_t
        %    tmp = spike_mat(1:num_neuron,t);
        %    tmp_fired = find(tmp>0);
        %    frate(t) = numel(tmp_fired);
        %    
        %    tmp = spike_mat(1:30:num_neuron,t);
        %    tmp_fired = find(tmp>0);   
        %    to_plot_y = [tmp_fired', to_plot_y];
        %    to_plot_x = [repmat(t,1,numel(tmp_fired)), to_plot_x];
        %    % scatter( repmat(t,1,numel(tmp_fired)), tmp_fired, '.' );
        %    hold on
        % end
        % scatter( to_plot_x, to_plot_y, '.' );
        % % p = pcolor(spike_mat(1:30:end,:)>0);
        % % set(p, 'EdgeColor', 'none')
        % xlim([0 num_draw_t]);
        % grid on
        % 
        % subplot(2,1,2)
        % plot(frate);
        % xlim([0 num_draw_t]);
        % grid on
        % 
        % 
        % 
        %% Frequency analysis
        % X = sum(spike_mat,1);
        % 
        % Fs = 1000;            % Sampling frequency
        % T = 1/Fs;             % Sampling period
        % L = size(X,2);        % Length of signal
        % t = (0:L-1)*T;        % Time vector
        % %plot(1000*t(1:50),X(1:50))
        % %title('Signal Corrupted with Zero-Mean Random Noise')
        % %xlabel('t (milliseconds)')
        % %ylabel('X(t)')
        % 
        % Y = fft(X);
        % P2 = abs(Y/L);
        % P1 = P2(1:L/2+1);
        % P1(2:end-1) = 2*P1(2:end-1);
        % f = Fs*(0:(L/2))/L;
        % 
        % figure
        % %plot(f,P1)
        % plot(f(2:end),P1(2:end))
        % xlabel('f (Hz)')
        % ylabel('|P1(f)|')
        % grid on    
    end
end

%%
f1 = figure;
f2 = figure;


fprintf('mean D2, for:\n')
fprintf('\n');
h = [];
name = {};
for id = 1:ntemplates
    
    figure(f1)
    m = mean(avg_max_d2(id,:), 2);
    s = std(avg_max_d2(id,:), [], 2);
    h(end+1) = errorbar(id, m, s, '.');
    hold on
    xlim([0 5])
    name{end+1} = strrep(dirnames_template{id}, 'network_full_leonardo_', '');
    name{end} = strrep(name{end}, 'Np_40_Ns_30_p_ratio_2', '');    
    name{end} = strrep(name{end}, '_rate100_run%d', '');
    name{end} = strrep(name{end}, '_edge_wrap_1_', '');
    name{end} = strrep(name{end}, '_', ' ');
    name{end} = strrep(name{end}, '/', ' ');
    
    fprintf('\t%s : %2.5f (%2.5f)\n', name{end}, m, s);
    
    figure(f2)
    subplot(2,2,id)
    histogram(avg_max_resolution(id,:), 5);
    title(name{end})    
end
fprintf('\n');
legend(h, name)


%%
f1 = figure;
f2 = figure;

fprintf('mean maximumal resolution, for:\n')
fprintf('\n');

h = [];
name = {};
for id = 1:ntemplates
    
    figure(f1)
    m = mean(avg_max_resolution(id,:), 2);
    s = std(avg_max_resolution(id,:), [], 2);
    h(end+1) = errorbar(id, m, s, '.');
    hold on
    xlim([0 5])
    ylim([min(avg_max_resolution(:))-100, max(avg_max_resolution(:)+100)])
    name{end+1} = strrep(dirnames_template{id}, 'network_full_leonardo_', '');
    name{end} = strrep(name{end}, 'Np_40_Ns_30_p_ratio_2', '');    
    name{end} = strrep(name{end}, '_rate100_run%d', '');
    name{end} = strrep(name{end}, '_edge_wrap_1_', '');
    name{end} = strrep(name{end}, '_', ' ');
    name{end} = strrep(name{end}, '/', ' ');
    
    fprintf('\t%s : %2.5f (%2.5f)\n', name{end}, m, s);
    
    figure(f2)
    subplot(2,2,id)
    histogram(avg_max_resolution(id,:), 5);
    title(name{end})
end

fprintf('\n');
legend(h, name)



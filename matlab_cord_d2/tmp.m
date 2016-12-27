

num_draw_t = 200;
num_draw_neuron = 1600;

frate = zeros(1,num_draw_t);
subplot(2,1,1)
for t=1:1:num_draw_t
   tmp = spike_mat(1:num_draw_neuron,t);
   tmp_fired = find(tmp>0);
   frate(t) = numel(tmp_fired); 
   scatter( repmat(t,1,numel(tmp_fired)), tmp_fired, '.' );
   hold on
end
xlim([0 num_draw_t]);
grid on

subplot(2,1,2)
plot(frate);
xlim([0 num_draw_t]);
grid on
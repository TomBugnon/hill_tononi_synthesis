dirname = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/examine_EPSP_IPSP_pablo/';

filename_Vp_v = strcat( dirname, 'connection_Vp_v_AMPA_syn.dat');
filename_Vp_h = strcat( dirname, 'connection_Vp_h_AMPA_syn.dat');
filename_Tp = strcat( dirname, 'connection_Tp_AMPA_syn.dat');

Vp_v = load(filename_Vp_v);
Vp_h = load(filename_Vp_h);
Tp = load(filename_Tp);

Vp_v_start_idx = min(Vp_v(:,1));
target_idx = Vp_v_start_idx + (3200+1600) + 800;

Vp_v_pre_rows = find(Vp_v(:,2)==target_idx);
Vp_h_pre_rows = find(Vp_h(:,2)==target_idx);
Tp_pre_rows = find(Tp(:,2)==target_idx);

Vp_v_pre = Vp_v(Vp_v_pre_rows,1);
Vp_h_pre = Vp_v(Vp_h_pre_rows,1);
Tp_pre = Vp_v(Tp_pre_rows,1);
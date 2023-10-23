%% 2d touch maap testing 5/23/23--> touch maps for S1/S2 for figure 2.
all_dat= get_s1s2_all_dat;
%vali= find(all_dat.has_effect_radius_data==1);
%cut= [1, 3:length(vali)];
%vali= vali([2:6 8]);
%dates_used= all_dat.lesion_dates; 
% S1 example
ani_str= 'an017510';
day_idx= 5;
area_nrni=  get_s1s2_neuron_subset_idx(all_dat.anim_data{23}.ids, 's1_all', all_dat, 23);
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)<10030000); %sub 1/6
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>10030000 & all_dat.anim_data{23}.ids(area_nrni)<10040000);%sub2
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>10040000 & all_dat.anim_data{23}.ids(area_nrni)<20000000);%sub3
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>20000000 & all_dat.anim_data{23}.ids(area_nrni)<20030000);    %sub 4
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>20030000 & all_dat.anim_data{23}.ids(area_nrni)<20040000);    %sub 5
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>20040000) ;    %sub 6

ai=23;
area_nrni_usw= intersect(all_dat.types_by_idx{ai}.usw_by_day{5}, area_nrni);
area_nrni_bsw=  intersect(all_dat.types_by_idx{ai}.bsw_by_day{5}, area_nrni);
area_nrni_mw= intersect(all_dat.types_by_idx{ai}.mw_by_day{5}, area_nrni);

area_ids_usw= area_nrni_usw(all_dat.anim_data{23}.ids(area_nrni_usw)<20030000);    %sub 4-6
area_ids_bsw= area_nrni_bsw(all_dat.anim_data{23}.ids(area_nrni_bsw)<20030000);    %sub 4-6
area_ids_mw= area_nrni_mw(all_dat.anim_data{23}.ids(area_nrni_mw)<20030000);    %sub 4-6
%area_ids= area_nrni_usw(all_dat.anim_data{23}.ids(area_nrni_usw)>10020000 & all_dat.anim_data{23}.ids(area_nrni_usw)<20000000);    %sub 1-3


fh=figure('Position',[0 0 600 600]);
    

    s1_w1p_ax = axes('Position', [0.1 0.55 0.4 0.4]);
    s2_w1p_ax = axes('Position', [0.55 0.55 0.4 0.4]);
    s1_w2p_ax = axes('Position', [0.1 0.1 0.4 0.4]);
    s2_w2p_ax = axes('Position', [0.55 0.1 0.4 0.4]);
    

    
    max_presp = 0.5;
    ball_size_range = [1 12];
    
    
    plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'usw', 'w1p',area_nrni_usw, 0, max_presp, ball_size_range); 
    plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'usw', 'w2p', area_nrni_usw, 0, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'bsw', 'w1p', area_nrni_bsw, 0, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'bsw', 'w2p', area_nrni_bsw, 0, max_presp, ball_size_range);

%     
   plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'mw', 'w1p', area_nrni_mw, 0, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'mw', 'w2p', area_nrni_mw, 0, max_presp, ball_size_range);
%     plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w2p_ax, 'mw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
%     
     %plot_2d_s1s2(all_dat, ani_str, day_post, s2_w1p_ax, 'usw', 'w1p', area_nrni, 1, max_presp, ball_size_range);
     %plot_2d_s1s2(all_dat, ani_str, day_post, s2_w1p_ax, 'bsw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
     %plot_2d_s1s2(all_dat, ani_str, day_post, s1_w1p_ax, 'usw', 'w1r', area_nrni, 0, max_presp, ball_size_range);
     %plot_2d_s1s2(all_dat, ani_str, day_post, s1_w1p_ax, 'usw', 'w2r', area_nrni, 0, max_presp, ball_size_range);

%     plot_2d_s1s2(all_dat, ani_str, day_post, s2_w1p_ax, 'mw', 'w1p', area_nrni, 0, max_presp, ball_size_range);
%     
   % plot_2d_s1s2(all_dat, ani_str, day_post, s2_w2p_ax, 'usw', 'w2p', area_nrni, 1, max_presp, ball_size_range);
%     plot_2d_s1s2(all_dat, ani_str, day_post, s2_w2p_ax, 'bsw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
%     plot_2d_s1s2(all_dat, ani_str, day_post, s2_w2p_ax, 'mw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
%     
%     
%     
        
%% s2
ani_str= 'an016650';
day_idx= 6;
area_nrni=  get_s1s2_neuron_subset_idx(all_dat.anim_data{20}.ids, 's2_all', all_dat, 20);
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)<10030000); %sub 1/6
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>10030000 & all_dat.anim_data{23}.ids(area_nrni)<10040000);%sub2
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>10040000 & all_dat.anim_data{23}.ids(area_nrni)<20000000);%sub3
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>20000000 & all_dat.anim_data{23}.ids(area_nrni)<20030000);    %sub 4
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>20030000 & all_dat.anim_data{23}.ids(area_nrni)<20040000);    %sub 5
%area_ids= area_nrni(all_dat.anim_data{23}.ids(area_nrni)>20040000) ;    %sub 6

ai=20;
area_nrni_usw= intersect(all_dat.types_by_idx{ai}.usw_by_day{6}, area_nrni);
area_nrni_bsw=  intersect(all_dat.types_by_idx{ai}.bsw_by_day{6}, area_nrni);
area_nrni_mw= intersect(all_dat.types_by_idx{ai}.mw_by_day{6}, area_nrni);

area_ids_usw= area_nrni_usw(all_dat.anim_data{20}.ids(area_nrni_usw)<1020030000);    %sub 4-6
area_ids_bsw= area_nrni_bsw(all_dat.anim_data{20}.ids(area_nrni_bsw)<1020030000);    %sub 4-6
area_ids_mw= area_nrni_mw(all_dat.anim_data{20}.ids(area_nrni_mw)<1020030000);    %sub 4-6
%area_ids= area_nrni_usw(all_dat.anim_data{23}.ids(area_nrni_usw)>10020000 & all_dat.anim_data{23}.ids(area_nrni_usw)<20000000);    %sub 1-3


fh2=figure('Position',[0 0 600 600]);
    

    s1_w1p_ax = axes('Position', [0.1 0.55 0.4 0.4]);
    s2_w1p_ax = axes('Position', [0.55 0.55 0.4 0.4]);
    s1_w2p_ax = axes('Position', [0.1 0.1 0.4 0.4]);
    s2_w2p_ax = axes('Position', [0.55 0.1 0.4 0.4]);
    

    
    max_presp = 0.5;
    ball_size_range = [1 12];
    
    
    plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'usw', 'w1p',area_nrni_usw, 0, max_presp, ball_size_range); 
    plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'usw', 'w2p', area_nrni_usw, 0, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'bsw', 'w1p', area_nrni_bsw, 0, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'bsw', 'w2p', area_nrni_bsw, 0, max_presp, ball_size_range);

%     
   plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'mw', 'w1p', area_nrni_mw, 0, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'mw', 'w2p', area_nrni_mw, 0, max_presp, ball_size_range); 
    
   
 
    
    
    %% 
   
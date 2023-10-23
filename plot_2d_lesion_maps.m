%% 2d touch maap testing 4/19/23
all_dat= get_s1s2_all_dat;
vali= find(all_dat.has_effect_radius_data==1);
%cut= [1, 3:length(vali)];
vali= vali([2:6 8]);
dates_used= all_dat.lesion_dates; 

for v= 1:length(vali)
    ai= vali(v); 
    ani_str = all_dat.anims{ai};
    dead_all= get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, 'rad_effect', all_dat, ai);
    area_nrni=  get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, 's1_manual_restrict', all_dat, ai);
    dead_id= intersect(dead_all, area_nrni); % dead id but only in manual restriction area 
    alive_id= setdiff(area_nrni, dead_id);
    
    pre_di = find(dates_used{ai} == -1);
    post_di = find(dates_used{ai} == 1);

        % map onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{ai}.date_str));
        for d=1:length(all_dat.anim_data{ai}.date_str)
            if (length(all_dat.anim_data{ai}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}));
            end
        end
        
        day_pre= find(d_num==pre_di(end));
        day_post= find(d_num==post_di(1)); 
        
      for i=2:2%1:3%3%%3 
        if i==1 area_ids= area_nrni(all_dat.anim_data{ai}.ids(area_nrni)<10030000); dead= dead_id(all_dat.anim_data{ai}.ids(dead_id)<10030000); alive=alive_id(all_dat.anim_data{ai}.ids(alive_id)<10030000); tp=all_dat.anim_data{ai}.ids(alive_id)<10030000; end 
        if i==2 area_ids= area_nrni(all_dat.anim_data{ai}.ids(area_nrni)>10030000 & all_dat.anim_data{ai}.ids(area_nrni)<10040000); dead= dead_id(all_dat.anim_data{ai}.ids(dead_id)>10030000 & all_dat.anim_data{ai}.ids(dead_id)<10040000); alive= alive_id(all_dat.anim_data{ai}.ids(alive_id)>10030000 & all_dat.anim_data{ai}.ids(alive_id)<10040000); tp= (all_dat.anim_data{ai}.ids(alive_id)>10030000 & all_dat.anim_data{ai}.ids(alive_id)<10040000); end 
        if i==3 area_ids= area_nrni(all_dat.anim_data{ai}.ids(area_nrni)<10040000) ; dead= dead_id(all_dat.anim_data{ai}.ids(dead_id)>10040000); alive=alive_id(all_dat.anim_data{ai}.ids(alive_id)>10040000); tp=all_dat.anim_data{ai}.ids(alive_id)>10040000; end 
    
        figure;
        x= all_dat.anim_data{ai}.x_um(dead);
        y= all_dat.anim_data{ai}.y_um(dead);
        plot(x,y,'ko');hold on; 
        xx= all_dat.anim_data{ai}.x_um(alive);
        yy= all_dat.anim_data{ai}.y_um(alive); 
        plot (xx, yy, 'go');
        title(['ani is ', ani_str])
        
      end 
      
    fh=figure('Position',[0 0 600 600]);

    s1_w1p_ax = axes('Position', [0.1 0.55 0.4 0.4]);
    s2_w1p_ax = axes('Position', [0.55 0.55 0.4 0.4]);
    s1_w2p_ax = axes('Position', [0.1 0.1 0.4 0.4]);
    s2_w2p_ax = axes('Position', [0.55 0.1 0.4 0.4]);

    %all_dat = get_s1s2_all_dat;
    %all_dat=all_dat;
    % parameters for plot_3d
    %ani_str = all_dat.anims{ai};%'an018920' % F4'an018927'; %an016650, F2 %S1 example an: an017510 day 5. S2 example an: 16650 6
    
    max_presp = 0.5;
    ball_size_range = [1 14];
    
    
    plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w1p_ax, 'usw', 'w1p',area_nrni, 1, max_presp, ball_size_range);
    plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w1p_ax, 'usw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
    %plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w1p_ax, 'usw', 'w1r', area_nrni, 0, max_presp, ball_size_range);
    %plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w1p_ax, 'usw', 'w2r', area_nrni, 0, max_presp, ball_size_range);

%     
     %plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w2p_ax, 'usw', 'w2p', area_nrni, 1, max_presp, ball_size_range);
%     plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w2p_ax, 'bsw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
%     plot_2d_s1s2(all_dat, ani_str, day_pre, s1_w2p_ax, 'mw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
%     
     plot_2d_s1s2(all_dat, ani_str, day_post, s2_w1p_ax, 'usw', 'w1p', area_nrni, 1, max_presp, ball_size_range);
     plot_2d_s1s2(all_dat, ani_str, day_post, s2_w1p_ax, 'bsw', 'w2p', area_nrni, 0, max_presp, ball_size_range);
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
        

end 
    
   
 
    
    
    %% 
   
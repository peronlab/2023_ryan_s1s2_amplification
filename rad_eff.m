%LR 4/19/23: testing script for determining radius of effect of lesion 

all_dat = get_s1s2_all_dat;
vali= find(all_dat.has_effect_radius_data==1);
%cut= [1:5, 7:length(vali)];
vali= vali([3:end]); %these 2 animals don't actually have enough left in the imaging FOV after lesion
dates_used= all_dat.lesion_dates;

%% 1. find neurons in a plane that are alive, dead
figure; 
%vali= [2 4]
for v= 1:length(vali)
    
    %ai=4; %test anim 17722
    ai= vali(v)
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
    
        for ii=1:2 
        if ii==1 dat_dayi= day_pre; end 
        if ii==2 dat_dayi= day_post;end 
        
        resp_binary_mat = nan*zeros(4,length(all_dat.anim_data{ai}.ids));
        resp_binary_mat(1,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w1p_by_day{dat_dayi})) = 1;
        resp_binary_mat(2,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w1r_by_day{dat_dayi})) = 1;
        resp_binary_mat(3,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w2p_by_day{dat_dayi})) = 1;
        resp_binary_mat(4,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w1r_by_day{dat_dayi})) = 1;

        resp_dff_mat = nan*zeros(4,length(all_dat.anim_data{ai}.ids));
        resp_dff_mat(1,:) = all_dat.anim_mats{ai}.meanRespW1p(dat_dayi,:);
        resp_dff_mat(2,:) = all_dat.anim_mats{ai}.meanRespW1r(dat_dayi,:);
        resp_dff_mat(3,:) = all_dat.anim_mats{ai}.meanRespW2p(dat_dayi,:);
        resp_dff_mat(4,:) = all_dat.anim_mats{ai}.meanRespW2r(dat_dayi,:);
        mu_resp_per_cell = nanmean(resp_dff_mat.*resp_binary_mat);
        
        if ii==1 pre= mu_resp_per_cell; end 
        if ii==2 post= mu_resp_per_cell; end 
    end 
         %% 2. for each alive neuron, find closest dead relative 
%    
% need to do this by plane to properly incorporate z here
    for i=1:3%%3 
    if i==1 dead= dead_id(all_dat.anim_data{ai}.ids(dead_id)<10030000); alive=alive_id(all_dat.anim_data{ai}.ids(alive_id)<10030000); tp=all_dat.anim_data{ai}.ids(alive_id)<10030000; end 
    if i==2 dead= dead_id(all_dat.anim_data{ai}.ids(dead_id)>10030000 & all_dat.anim_data{ai}.ids(dead_id)<10040000); alive= alive_id(all_dat.anim_data{ai}.ids(alive_id)>10030000 & all_dat.anim_data{ai}.ids(alive_id)<10040000); tp= (all_dat.anim_data{ai}.ids(alive_id)>10030000 & all_dat.anim_data{ai}.ids(alive_id)<10040000); end 
    if i==3 dead= dead_id(all_dat.anim_data{ai}.ids(dead_id)>10040000); alive=alive_id(all_dat.anim_data{ai}.ids(alive_id)>10040000); tp=all_dat.anim_data{ai}.ids(alive_id)>10040000; end 
    
    ani_numalive (v,i)= length(alive); 
    
    pre_mu= pre(tp);
    post_mu= post(tp);
    
    %P(:,1)= all_dat.anim_data{ai}.x_um(dead); % what you're comparing to-- DeAD (making a matrix thats x, y) 
    %P(:,2)= all_dat.anim_data{ai}.y_um(dead);
    %PQ (:,1)= all_dat.anim_data{ai}.x_um(alive);
    %PQ (:,2)= all_dat.anim_data{ai}.y_um(alive); %what you are querying through -- ALIVE
    
    xd= all_dat.anim_data{ai}.x_um(dead);
    yd= all_dat.anim_data{ai}.y_um(dead);
    for nn= 1:length(xd)
        P(nn,:)= [xd(nn),yd(nn)];
    end 
    
    xa= all_dat.anim_data{ai}.x_um(alive);
    ya= all_dat.anim_data{ai}.y_um(alive);
    for mm=1:length(xa)
         PQ(mm,:)= [xa(mm),ya(mm)];
    end 
    [closest_dead, dist]= dsearchn(P,PQ); %this looks at x,y coord from each alive cell and finds the nearest dead
    b= 50;
    bins= [0:b:400]; 
    for bb=1: (length(bins)-1)
        idxes= find(dist> bins(bb) & dist < bins(bb+1));
        alive_thisbin= alive(idxes);
        preresp= pre_mu(idxes);
        postresp= post_mu(idxes);
        resp_this_bin_pre(bb,i)= sum(preresp(~isnan(preresp)));
        resp_this_bin_post(bb,i)= sum(postresp(~isnan(postresp)));
        
        num_touch_pre(bb,i)= length(intersect(all_dat.types_by_idx{ai}.usw_by_day{day_pre}, alive_thisbin))+ length(intersect(all_dat.types_by_idx{ai}.bsw_by_day{day_pre}, alive_thisbin))+ length(intersect(all_dat.types_by_idx{ai}.mw_by_day{day_pre}, alive_thisbin));
        num_touch_post(bb,i)= length(intersect(all_dat.types_by_idx{ai}.usw_by_day{day_post}, alive_thisbin))+ length(intersect(all_dat.types_by_idx{ai}.bsw_by_day{day_post}, alive_thisbin))+ length(intersect(all_dat.types_by_idx{ai}.mw_by_day{day_post}, alive_thisbin));

        
        num_touch_pre_W1P(bb,i)= length (intersect(all_dat.types_by_idx{ai}.touch_w1p_by_day{day_pre}, alive_thisbin))+length (intersect(all_dat.types_by_idx{ai}.touch_w1r_by_day{day_pre}, alive_thisbin));
        num_touch_post_W1P(bb,i)= length (intersect(all_dat.types_by_idx{ai}.touch_w1p_by_day{day_post}, alive_thisbin))+length (intersect(all_dat.types_by_idx{ai}.touch_w1r_by_day{day_post}, alive_thisbin));
        num_touch_pre_W2P (bb,i)=  length (intersect(all_dat.types_by_idx{ai}.touch_w2p_by_day{day_pre}, alive_thisbin))+length (intersect(all_dat.types_by_idx{ai}.touch_w2r_by_day{day_pre}, alive_thisbin));
        num_touch_post_W2P(bb,i)= length (intersect(all_dat.types_by_idx{ai}.touch_w2p_by_day{day_post}, alive_thisbin))+length (intersect(all_dat.types_by_idx{ai}.touch_w2r_by_day{day_post}, alive_thisbin));
        
         
        for c=1:2 
            if c==1 day_idx= day_pre; end 
            if c==2 day_idx= day_post; end 
            
            w1p_idx= (intersect(all_dat.types_by_idx{ai}.touch_w1p_by_day{day_idx}, alive_thisbin));
            w2p_idx= (intersect(all_dat.types_by_idx{ai}.touch_w2p_by_day{day_idx}, alive_thisbin));
            w1r_idx= (intersect(all_dat.types_by_idx{ai}.touch_w1r_by_day{day_idx}, alive_thisbin));
            w2r_idx= (intersect(all_dat.types_by_idx{ai}.touch_w2r_by_day{day_idx}, alive_thisbin));
            
            w1p_resp= (all_dat.anim_mats{ai}.meanRespW1p(day_idx,:));
            w2p_resp= (all_dat.anim_mats{ai}.meanRespW2p(day_idx,:));
            w1r_resp= (all_dat.anim_mats{ai}.meanRespW1r(day_idx,:));
            w2r_resp= (all_dat.anim_mats{ai}.meanRespW2r(day_idx,:));
            
            w1p(bb,i,c)= sum(w1p_resp(w1p_idx));
            w2p(bb,i,c)= sum(w2p_resp(w2p_idx)); 
            w1r(bb,i,c)= sum(w1r_resp(w1r_idx)); 
            w2r(bb,i,c)= sum(w1r_resp(w2r_idx));
            
            %if c==1 w1p_pre= w1p; w2p_pre=w2p; w1r_pre=w1r; w2r_pre=w2r; end 
            %if c==2 w1p_post= w1p; w2p_post= w2p; w1r_post= w1r; w2r_post= w2r; end 
            
            %clear w1p; clear w2p; clear w1r; clear w2r;
        end
        
           

%         num_cells_pre(1, bb,v)= length(intersect(all_dat.types_by_idx{ai}.usw_by_day{day_pre}, alive_thisbin));
%         num_cells_post(1, bb,v)= length(intersect(all_dat.types_by_idx{ai}.usw_by_day{day_post}, alive_thisbin));
%         num_cells_pre(2, bb,v)= length(intersect(all_dat.types_by_idx{ai}.bsw_by_day{day_pre}, alive_thisbin));
%         num_cells_post(2, bb,v)= length(intersect(all_dat.types_by_idx{ai}.bsw_by_day{day_post}, alive_thisbin));
%         num_cells_pre(3, bb,v)= length(intersect(all_dat.types_by_idx{ai}.mw_by_day{day_pre}, alive_thisbin));
%         num_cells_post(3, bb,v)= length(intersect(all_dat.types_by_idx{ai}.mw_by_day{day_post}, alive_thisbin));
%       
    end 
    
      
        
    
    clear P; clear PQ; 
   


    end
  
    pre_nt(v,:)= sum(num_touch_pre,2);
    post_nt(v,:)= sum(num_touch_post,2);
    
    sum_w1p= sum(w1p,2);
    sum_w2p= sum(w2p,2);
    sum_w1r= sum(w1r, 2);
    sum_w2r= sum(w2r,2);
    
    ani_per_resp_W1p (v,:)= sum_w1p(:,:,2)- sum_w1p(:,:,1);
    ani_per_resp_W1r (v,:)= sum_w1r(:,:,2)- sum_w1r(:,:,1); 
    ani_per_resp_W2p (v,:)= sum_w2p(:,:,2)- sum_w2p(:,:,1); 
    ani_per_resp_W2r (v,:)= sum_w2r(:,:,2)- sum_w2r(:,:,1); 

      
        if strcmp(ani_str, 'an017722') || strcmp(ani_str, 'an014355')
            lesioned_nt_pre(v,:)= sum(num_touch_pre_W1P,2); 
            lesioned_nt_post(v,:)= sum(num_touch_post_W1P,2); 
            unlesion_nt_pre (v,:)= sum(num_touch_pre_W2P,2);
            unlesion_nt_post (v,:)= sum(num_touch_post_W2P,2);
            lesioned_pro (v,:)= ani_per_resp_W1p(v,:);
            lesioned_ret (v,:)= ani_per_resp_W1r(v,:);
            unlesion_pro (v,:)= ani_per_resp_W2p(v,:);
            unlesion_ret(v,:)= ani_per_resp_W2r(v,:); 
        else 
            lesioned_nt_pre(v,:)= sum(num_touch_pre_W2P,2); 
            lesioned_nt_post(v,:)= sum(num_touch_post_W2P,2); 
            unlesion_nt_pre (v,:)= sum(num_touch_pre_W1P,2);
            unlesion_nt_post (v,:)= sum(num_touch_post_W1P,2);
            lesioned_pro (v,:)= ani_per_resp_W2p (v,:);
            lesioned_ret (v,:)= ani_per_resp_W2r (v,:);
            unlesion_pro (v,:)= ani_per_resp_W1p (v,:);
            unlesion_ret(v,:)= ani_per_resp_W1r (v,:); 
        end
    
     
     
   

   
end 
%% plotting 

    x_ax= bins(2:end);%[25 75 125 175 225 275 325 375];
    color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
   
    for bh= 1:length(vali)
        
        new= (post_nt(bh,:));
        old= (pre_nt(bh,:));
        diffy(bh,:)= new -old; 
        new_les= lesioned_nt_post(bh,:);
        old_les= lesioned_nt_pre(bh,:);
        diffy_les (bh,:)= new_les-old_les;
        new_unles= unlesion_nt_post(bh,:);
        old_unles=  unlesion_nt_pre(bh,:);
        diffy_unles(bh,:)= new_unles-old_unles;
        
        plot (x_ax-(b/2), new_les-old_les, 'b-'); hold on;
        plot (x_ax-(b/2), new_unles-old_unles, 'r-');hold on;
       
    end 
    
    plot (x_ax-(b/2), nanmean(diffy), 'k-', 'LineWidth', 7); hold on;
    plot (x_ax-(b/2), nanmean(diffy_les), 'b-', 'LineWidth', 5); 
    plot (x_ax-(b/2), nanmean(diffy_unles), 'r-', 'LineWidth', 5); 
    xlim([0 400])
    yline(0,'k--')
    
    maxi_l= max(abs(lesioned_pro)');
    maxi_u= max(abs(unlesion_pro)');
    range_l= range(lesioned_pro,2);
    range_u= range(unlesion_pro,2);
    range_perani= (range_l+range_u)/2;
   
    figure (2);
     for bbh= 1:length(vali)
          normd_lp (bbh,:)= lesioned_pro(bbh,:)/range_perani(bbh);
          normd_up (bbh,:)= unlesion_pro(bbh,:)/range_perani(bbh);
         
          plot (x_ax-(b/2), lesioned_pro(bbh,:)/range_perani(bbh),'r-'); hold on;
          plot(x_ax-(b/2), unlesion_pro(bbh,:)/range_perani(bbh), 'k-'); 
        % plot (x_ax-(50/2), lesioned_pro(bbh,:),'r-'); hold on;
         %plot(x_ax-(50/2), unlesion_pro(bbh,:), 'k-'); 
     end 
     %ylim([-12 5])
     plot (x_ax-(b/2), nanmean(normd_lp),'r-', 'LineWidth',5); hold on;
     %plot (x_ax-(50/2), nanmean(normd_lp),'r-', 'LineWidth',5); hold on;
     %plot (x_ax-(50/2), nanmean(normd_up),'k-', 'LineWidth',5); 
   %plot (x_ax-(50/2), nanmean(lesioned_ret), 'r--', 'LineWidth',5); 
    plot (x_ax-(b/2), nanmean(normd_up), 'k-', 'LineWidth',5); 
    yline(0,'k--')
     %plot (x_ax- (50/2), nanmean(unlesion_ret),'k--', 'LineWidth', 5);
  
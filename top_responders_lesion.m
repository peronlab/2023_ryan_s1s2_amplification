%% LR new script to test if responsiveness of cells changes pre/post lesion. 
    
    all_dat = get_s1s2_all_dat;
    s1_lesion_ani = find(~strcmp('',all_dat.s1_lesion_date));
    s2_lesion_ani = find(~strcmp('',all_dat.s2_lesion_date));
%   sham_lesion_ani = find(all_dat.is_1daysham); % restrict to only 1day shams (This is more appropriate - what we do for everything else!)
    sham_lesion_ani = find(~strcmp('',all_dat.sham_lesion_date)); % use all shams
    
    fh = figure ('Position', [0 400 600 400]);
    ax1 = subplot('Position', [.1 .1 .4 .4]);hold on; 
    ax2= subplot('Position', [.1 .55 .4 .4]);hold on;
    
    ax3= subplot('Position', [.6 .55 .4 .4]);hold on;
    ax4= subplot('Position', [.6 .1 .4 .4]);hold on;
    %things to put in the function when we make it 
    %% 
    plot_resp_per_cell_prepost_lesion(all_dat,1, ax1, s1_lesion_ani, 's2_manual_restrict', all_dat.lesion_dates, 'S2 after S1 lesion, neurons that keep ID'); 
    plot_resp_per_cell_prepost_lesion(all_dat,1, ax2, s2_lesion_ani, 's1_manual_restrict', all_dat.lesion_dates, 'S1 after S2 lesion, neurons that keep ID'); 
    %plot_resp_per_cell_prepost_lesion(all_dat, 1, ax3, sham_lesion_ani, 's1_manual restrict', all_dat.sham_dates, 'S1 after sham');
    [num_cellsham, num_turnsham, per_tsham, per_bsham]= plot_top_v_bottom(all_dat, ax3, sham_lesion_ani, 's1_manual_restrict', all_dat.sham_dates, 'S1 after saahm  lesion, top 10per v bottom 50% ID keep'); 
    [num_cells, num_turn, per_tS2, per_bS2]= plot_top_v_bottom(all_dat, ax4, s2_lesion_ani, 's1_manual_restrict', all_dat.lesion_dates, 'S1 after S2per_ lesion, top 10per v bottom 50% ID keep'); 
    
    %plot_resp_per_cell_prepost_lesion(all_dat,2, ax3, s1_lesion_ani, 's2_manual_restrict', all_dat.lesion_dates, 'S2 after S1 lesion, ID Change'); 

%%
%% goal: 1. find top 10% of MW cells. 2. how many of those keep their identity, how many shift? bottom 50%.

function [num_cells, num_turnover, per_top_stay, per_bottom_stay]= plot_top_v_bottom (all_dat, ax, vali, restrict_tag, dates_used, tstr)
 for v=1:length(vali)
    ai=vali(v);
    area_nrni= get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, restrict_tag, all_dat, ai); 
    num_cells(v)= length(area_nrni);
    pre_di = find(dates_used{ai} == -1);
    post_di = find(dates_used{ai} == 1);

        % map onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{ai}.date_str));
        for d=1:length(all_dat.anim_data{ai}.date_str)
            if (length(all_dat.anim_data{ai}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}));
            end
        end
        
    
    
    for i= 1:2 %% make a matrix for pre, then post, values 
        if i==1
            for di= 1:length(pre_di)
                dat_dayi = find(d_num == pre_di(di));
            end 
       end 
        if i==2  
           for di= 1:length(post_di)
                dat_dayi = find(d_num == post_di(di));
            end 
        end 
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
            mu_resp_per_cell = nanmean((resp_dff_mat.*resp_binary_mat));
            
            mw_idx= intersect(area_nrni,all_dat.types_by_idx{ai}.mw_by_day{dat_dayi}); 
            bsw_idx= intersect(area_nrni, all_dat.types_by_idx{ai}.bsw_by_day{dat_dayi});
            usw_idx= intersect(area_nrni, all_dat.types_by_idx{ai}.usw_by_day{dat_dayi});
            
            if i==1 
                mw_pre= mw_idx; 
                bsw_pre= bsw_idx; 
                usw_pre= usw_idx; 
                usw_resp_pre =nanmean(resp_dff_mat); 
                pre_resp= mu_resp_per_cell; 
            end
            
            if i==2
                mw_post= mw_idx; 
                bsw_post= bsw_idx; 
                usw_post= usw_idx; 
                usw_resp_post = nanmean(resp_dff_mat); 
                post_resp= mu_resp_per_cell; 
            end  
           
    end 
 
 
 for cti=1:3 
     if cti==1 pre=usw_pre; post= usw_post; resp= usw_resp_pre; end 
     if cti==2 pre= bsw_pre; post= bsw_post; resp= pre_resp; end 
     if cti==3 pre= mw_pre; post= mw_post; resp= pre_resp; end 
     
    top_per= 0.1; %starting by looking at top 10%
    bottom_per= 0.50;%bottom 50% 
    cells_vals= resp(pre); 
    if length(cells_vals)>1 
        sort_me= sort(cells_vals,'descend'); 
        top_resp= sort_me(1:round(length(cells_vals)*top_per)); %these are top x% of means of MW cells 
        bot_resp= sort_me(end-(round(length(cells_vals)*bottom_per)):end);
        idx_top= find(ismember(cells_vals, top_resp)); 
        idx_bottom= find(ismember(cells_vals,bot_resp));
        idx_stay= find(ismember(pre,post)); 
        per_top_stay(v, cti)= (length(intersect(idx_top,idx_stay))/length(idx_top)); 
        per_bottom_stay(v,cti)= (length(intersect(idx_bottom,idx_stay))/length(idx_bottom)); 
        num_turnover (v,cti)= (length(pre)- length(idx_stay))/length(pre);
    else 
       if ismember(pre,post)
           per_top_stay(v,cti)= 1; 
           per_bottom_stay(v,cti)= 0; 
       else 
           per_top_stay(v,cti)= nan; 
           per_bottom_stay(v, cti)= nan; 
       end 
    end 
 end 
 end 
 
 color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
 x= [1 3 5]; 

 x2= [2 4 6]; 
for i=1:3 
    bar(ax, x(i), nanmean(per_top_stay(:,i)),'FaceColor', color_arr(i,:)); hold on; 
    plot (ax, x(i), per_top_stay(:,i), 'ko');
    bar (ax, x2(i), nanmean(per_bottom_stay(:,i)), 'FaceColor', color_arr(i,:));
    plot(ax, x2(i), per_bottom_stay(:,i), 'ko');
    %bar(ax, x(i), nanmean(num_turnover(:,i)), 'FaceColor',
    %color_arr(i,:)); hold on; %this is for turnover rate supp 
    title (ax,tstr)
    ax.YLim =([0 1.1])
      
end 
end 
%% 

function plot_resp_per_cell_prepost_lesion(all_dat,type, ax, vali, restrict_tag, dates_used, tstr)

for v=1:length(vali)
    ai=vali(v);
    area_nrni= get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, restrict_tag, all_dat, ai); 
    pre_di = find(dates_used{ai} == -1);
    post_di = find(dates_used{ai} == 1);

        % map onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{ai}.date_str));
        for d=1:length(all_dat.anim_data{ai}.date_str)
            if (length(all_dat.anim_data{ai}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}));
            end
        end
        
    %day_pre= find(dates_used{ai} == -1);
    %day_pre= day_pre(end); %for now lets do day before, day after only 
    %day_post = find(dates_used{ai} == 1);
    %day_post= day_post(1);  
    
    for i= 1:2 %% make a matrix for pre, then post, values 
        if i==1
            for di= 1:length(pre_di)
                dat_dayi = find(d_num == pre_di(di));
            end 
       end 
        if i==2  
           for di= 1:length(post_di)
                dat_dayi = find(d_num == post_di(di));
            end 
        end 
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
            mu_resp_per_cell = nanmean((resp_dff_mat.*resp_binary_mat));
            
            mw_idx= intersect(area_nrni,all_dat.types_by_idx{ai}.mw_by_day{dat_dayi}); 
            bsw_idx= intersect(area_nrni, all_dat.types_by_idx{ai}.bsw_by_day{dat_dayi});
            usw_idx= intersect(area_nrni, all_dat.types_by_idx{ai}.usw_by_day{dat_dayi});
            
            if i==1 
                mw_pre= mw_idx; 
                bsw_pre= bsw_idx; 
                usw_pre= usw_idx; 
                usw_resp_pre =nanmean(resp_dff_mat); 
                pre_resp= mu_resp_per_cell; 
            end
            
            if i==2
                mw_post= mw_idx; 
                bsw_post= bsw_idx; 
                usw_post= usw_idx; 
                usw_resp_post = nanmean(resp_dff_mat); 
                post_resp= mu_resp_per_cell; 
            end  
           
    end 
    
    if type==1
     mw_cells_static(:,1)= pre_resp(intersect(mw_pre,mw_post));
     mw_cells_static (:,2)= post_resp(intersect(mw_pre, mw_post));
     bsw_cells_static (:,1)= pre_resp(intersect(bsw_pre, bsw_post));
     bsw_cells_static (:,2)= post_resp(intersect(bsw_pre, bsw_post));
     usw_cells_static (:,1)= usw_resp_pre(intersect(usw_pre, usw_post)); 
     usw_cells_static (:,2)= usw_resp_post(intersect(usw_pre, usw_post)); 
    end 
    % mw_response_pre= pre_resp(mw_pre); mw_response(:,2)= post_resp(mw_pre);
     
    color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];

    if type==1
        plot (ax, mw_cells_static(:,1), mw_cells_static(:,2), 'o', 'Color', color_arr(3,:)); hold on; 
        plot(ax, [-0.2 1.5],[-0.2 1.5],'--k');      
        plot (ax, bsw_cells_static (:,1), bsw_cells_static(:,2), 'o', 'Color', color_arr(2,:))
        plot (ax, usw_cells_static (:,1), usw_cells_static (:,2), 'o', 'Color', color_arr(1,:))
        ax.YLim =([0 1.5])
        ax.XLim =([0 1.5])
        hold on; 
        title (ax,tstr)
        clear mw_cells_static; clear bsw_cells_static; clear usw_cells_static;
    end 
    
    
    
end
end 
    
    
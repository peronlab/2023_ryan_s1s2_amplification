%
% Figure for following panels:
%   - fraction of all cells by area that have specific type (uSW, bSW, MW)
%   - mean encoding score for touch cell (may want to take top 5% and mean of that, for example)
%   - ??? 
% FOR NOW: superficial and deep as 2 seperate figs. removed 17518 b/c touch
% cell count super low. removed 14359 b/c weird in terms of
% ssubvolumes/nonmerging. trying to figure out a better way to handle. N= 7
% ffor now. 

function figure_compare_s1_s2_touch
    figure_compare_s1_s2_touch_super;
    figure_compare_s1_s2_touch_deep;
end 

function figure_compare_s1_s2_touch_super
    
    all_dat = get_s1s2_all_dat;
    anims_used = {'an014351','an016623','an017517','an016650','an016652','an014363','an017510'};
    ani_used = find(ismember(all_dat.anims, anims_used));
    %ani_used = [1  10    12    18    19    20    21]%17518 is weird, maybe remove. so is 14359
%find(all_dat.has_deep_data_s1s2);

    % --- make the blank plots
    fh = figure('Position', [0 0 2000 500]); 
    
    frac_by_type_and_area_ax = subplot('Position', [.05 .05 .15 .8]);  
    hold on;
    ax_2 = subplot('Position', [.25 .05 .15 .8]);  
    hold on;
    ax_mean= subplot('Position', [.45 .05 .15 .8]);   
    hold on;
    frac_all_ax = subplot('Position', [.65 .05 .15 .8]);  
    hold on;
   
    
    % --- fraction of cells by area and type
    frac_all_cells_by_type_and_area = nan*zeros(length(all_dat.anims), 2, 3); % animal X s1/s2 X uSW/bSW/MW
    frac_touch_cells_by_type_and_area = nan*zeros(length(all_dat.anims), 2, 3); % animal X s1/s2 X uSW/bSW/MW
    mean_S1_to_type= nan*zeros(4,length(all_dat.anims),3);
    mean_S2_to_type= nan*zeros(4,length(all_dat.anims),3);
    prob_S1_to_type= nan*zeros(4,length(all_dat.anims),3);
    prob_S2_to_type= nan*zeros(4,length(all_dat.anims),3);
    %mean_prob_resp_w1p= nan*zeros(length(all_dat.anims),2,3); %animal, s1/s2, usw/bsw/mw -- this will contain mean probability of response to w1p for each cell type in s1 and s2
    %S1_MW_resp =nan*zeros(100,4,length(all_dat.anims))
    
    d_s1 = 1; % most superficial of merged subvolumes
    d_s2 = 3; 
    d_s1_deep=2;
    d_s2_deep=4; 
    
    
    for ai=1:length(ani_used)
        a = ani_used(ai);

        % figure out which neurons go where
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);
        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
%        s1i = find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s1_id_range{a}(2));
%        s2i = find(all_dat.anim_data{a}.ids > all_dat.s2_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s2_id_range{a}(2));

        % which date should we use? second non-compound date
        non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);

        d_s1_used = d_s1;
        d_s2_used = d_s2;
       
        % animal-specific stuff
        if (strcmp(all_dat.anims{a}, 'an014359')) %has no merging, so need to use the best day, but also need to seperate superficial and deep!
            d_s2_used = 1;
        
%            s1i= s1i(all_dat.anim_data{3}.ids(s1i)<20000000);

            %s1i_deep= s1i(all_dat.anim_data{3}.ids(s1i)>20000000);
%            s2i= s2i(all_dat.anim_data{3}.ids(s2i)<920000000);
            %s2i_deep= s2i(all_dat.anim_data{3}.ids(s2i)>920000000);
        elseif (strcmp(all_dat.anims{a}, 'an016652'))
            d_s1_used = 3;
            d_s2_used = 1;
            
        end
        
       
        frac_all_cells_by_type_and_area (a, 1, 1) = length(intersect(s1i,all_dat.types_by_idx{a}.usw_by_day{d_s1_used}))/length(s1i);
        frac_all_cells_by_type_and_area (a, 2, 1) = length(intersect(s2i,all_dat.types_by_idx{a}.usw_by_day{d_s2_used}))/length(s2i);

        frac_all_cells_by_type_and_area (a, 1, 2) = length(intersect(s1i,all_dat.types_by_idx{a}.bsw_by_day{d_s1_used}))/length(s1i);
        frac_all_cells_by_type_and_area (a, 2, 2) = length(intersect(s2i,all_dat.types_by_idx{a}.bsw_by_day{d_s2_used}))/length(s2i);

        frac_all_cells_by_type_and_area (a, 1, 3) = length(intersect(s1i,all_dat.types_by_idx{a}.mw_by_day{d_s1_used}))/length(s1i);
        frac_all_cells_by_type_and_area (a, 2, 3) = length(intersect(s2i,all_dat.types_by_idx{a}.mw_by_day{d_s2_used}))/length(s2i);

        all_touch_s1_idx = [all_dat.types_by_idx{a}.usw_by_day{d_s1_used} all_dat.types_by_idx{a}.bsw_by_day{d_s1_used} all_dat.types_by_idx{a}.mw_by_day{d_s1_used}];
        all_touch_s2_idx = [all_dat.types_by_idx{a}.usw_by_day{d_s2_used} all_dat.types_by_idx{a}.bsw_by_day{d_s2_used} all_dat.types_by_idx{a}.mw_by_day{d_s2_used}];
       
        frac_touch_cells_by_type_and_area (a, 1, 1) = length(intersect(s1i,all_dat.types_by_idx{a}.usw_by_day{d_s1_used}))/length(all_touch_s1_idx);
        frac_touch_cells_by_type_and_area (a, 2, 1) = length(intersect(s2i,all_dat.types_by_idx{a}.usw_by_day{d_s2_used}))/length(all_touch_s2_idx);

        frac_touch_cells_by_type_and_area (a, 1, 2) = length(intersect(s1i,all_dat.types_by_idx{a}.bsw_by_day{d_s1_used}))/length(all_touch_s1_idx);
        frac_touch_cells_by_type_and_area (a, 2, 2) = length(intersect(s2i,all_dat.types_by_idx{a}.bsw_by_day{d_s2_used}))/length(all_touch_s2_idx);

        frac_touch_cells_by_type_and_area (a, 1, 3) = length(intersect(s1i,all_dat.types_by_idx{a}.mw_by_day{d_s1_used}))/length(all_touch_s1_idx);
        frac_touch_cells_by_type_and_area (a, 2, 3) = length(intersect(s2i,all_dat.types_by_idx{a}.mw_by_day{d_s2_used}))/length(all_touch_s2_idx);
            
        % make matrix of MW Responsive cells
            mw_idx_S1= intersect(all_dat.types_by_idx{a}.mw_by_day{d_s1_used},s1i); %find ids of multiwhisker cells in each subvolume
            mw_idx_S2= intersect(all_dat.types_by_idx{a}.mw_by_day{d_s2_used},s2i);
            bw_idx_S1= intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s1_used},s1i);
            bw_idx_S2= intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s2_used},s2i);

            w1p_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d_s1_used}, s1i); %find indexes of each type of touch in s1 and s2
            w1r_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d_s1_used}, s1i);
            w2p_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d_s1_used}, s1i);
            w2r_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d_s1_used}, s1i);    

            w1p_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d_s2_used}, s2i);
            w1r_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d_s2_used}, s2i);
            w2p_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d_s2_used}, s2i);
            w2r_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d_s2_used}, s2i);  

            for i= 1: length(s1i) %loop thorugh the possible all cell ids and make a matrix of cell id x type (w1p, w1r, w2p, w2r) per animal in S1 
               touch_types_within_S1(i,1, a)= any(w1p_idx_S1(:)== s1i(i));
               touch_types_within_S1(i,2,a)= any(w1r_idx_S1(:)== s1i(i));
               touch_types_within_S1(i,3,a)= any(w2p_idx_S1(:)== s1i(i));
               touch_types_within_S1(i, 4,a)= any(w2r_idx_S1(:)== s1i(i));   
            end      

             for i= 1: length(s2i) % same for S2
               touch_types_within_S2(i,1, a)= any(w1p_idx_S2(:)== s2i(i));
               touch_types_within_S2(i,2,a)= any(w1r_idx_S2(:)== s2i(i));
               touch_types_within_S2(i,3,a)= any(w2p_idx_S2(:)== s2i(i));
               touch_types_within_S2(i, 4,a)= any(w2r_idx_S2(:)== s2i(i));   
             end 
            
        % so touch types within S1 gives for each cell whether it was
        % response to w1p, then w1r, then w2p, then w2r-- full matrix, can
        % subselect by subtype of touch later. 
        
        S1_resp_mean= nan*zeros(max(s1i),4);
        S2_resp_mean= nan*zeros(max(s2i),4);
        S1_resp_prob= nan*zeros(max(s1i),4);
        S2_resp_prob= nan*zeros(max(s2i),4);
        %responses to each touch type (mean)
        for tt= 1:4 % loop through 4 touch types: w1p, w1r, w2p, w2r
           if tt==1 response=all_dat.types_by_idx{a}.ev_touch_w1p_by_day, mean_use=all_dat.anim_mats{a}.meanRespW1p; prob_use=all_dat.anim_mats{a}.probRespW1p; end 
           if tt==3 response=all_dat.types_by_idx{a}.ev_touch_w2p_by_day; mean_use= all_dat.anim_mats{a}.meanRespW2p; prob_use=all_dat.anim_mats{a}.probRespW2p;end 
           if tt==2 response=all_dat.types_by_idx{a}.ev_touch_w1r_by_day; mean_use= all_dat.anim_mats{a}.meanRespW1r;prob_use=all_dat.anim_mats{a}.probRespW1r;end 
           if tt==4 response=all_dat.types_by_idx{a}.ev_touch_w2r_by_day; mean_use= all_dat.anim_mats{a}.meanRespW2r;prob_use=all_dat.anim_mats{a}.probRespW2r;end 
          
           %example: mean_S1_to_type(1,:,1) is
           %the mean response in S1 for all animals to of usw cells to w1p.
           %this works for USw but not hte others, need to do those
           %seperately. 
           %order: usw,bw, mw
           mean_S1_to_type(tt,a,1)= mean(mean_use(d_s1_used, intersect(response{d_s1_used},all_dat.types_by_idx{a}.usw_by_day{d_s1_used})));
           
           mean_S2_to_type(tt, a, 1)= mean(mean_use(d_s2_used, intersect(response{d_s2_used},all_dat.types_by_idx{a}.usw_by_day{d_s2_used})));
           
           prob_S1_to_type(tt,a,1)= mean(prob_use(d_s1_used, intersect(response{d_s1_used},all_dat.types_by_idx{a}.usw_by_day{d_s1_used})));
           
           prob_S2_to_type(tt, a, 1)= mean(prob_use(d_s2_used, intersect(response{d_s2_used},all_dat.types_by_idx{a}.usw_by_day{d_s2_used})));
           
          % for each animal, construct a matrix of mean response of each
          % neuron to each touch type. we can take the dot product of this
          % matrix with the matrix listing what is resposive to what 
         
          S1_resp_mean (1:max(s1i),tt)= mean_use(d_s1_used, s1i);
          S2_resp_mean (min(s2i):max(s2i),tt)= mean_use(d_s2_used, s2i);
          S1_resp_prob (1:max(s1i),tt)= prob_use(d_s1_used, s1i);
          S2_resp_prob (min(s2i):max(s2i),tt)= prob_use(d_s2_used, s2i);
          
          S1_resp_mean_MW= S1_resp_mean(mw_idx_S1,:); %use only MW cells for dot prod.
          S1_resp_prob_MW= S1_resp_prob(mw_idx_S1,:); %use only MW cells for dot prod.

          touch_types_S1_MW= touch_types_within_S1(mw_idx_S1, :,a); %use only MW cells for dot prod
          
          S2_resp_mean_MW= S2_resp_mean(mw_idx_S2,:); %
          S2_resp_prob_MW= S2_resp_prob(mw_idx_S2,:); %

          touch_types_S2_MW= touch_types_within_S2(mw_idx_S2-max(s1i),:,a);
          
         
          mean_each_type_MW (a,:,1)= dot(S1_resp_mean_MW, touch_types_S1_MW)/length(touch_types_S1_MW); %S1
          mean_each_type_MW(a,:,2)= dot(S2_resp_mean_MW, touch_types_S2_MW)/length(touch_types_S2_MW); %S2
          
          prob_each_type_MW (a,:,1)= dot(S1_resp_prob_MW, touch_types_S1_MW)/length(touch_types_S1_MW); %S1
          prob_each_type_MW(a,:,2)= dot(S2_resp_prob_MW, touch_types_S2_MW)/length(touch_types_S2_MW); %S2
          
          %ok now do same dot product thing for bidir cells
          S1_resp_mean_bW= S1_resp_mean(bw_idx_S1,:); %use only MW cells for dot prod.
          S1_resp_prob_bW= S1_resp_prob(bw_idx_S1,:);
          
          touch_types_S1_bW= touch_types_within_S1(bw_idx_S1, :,a); %use only MW cells for dot prod
          S2_resp_mean_bW= S2_resp_mean(bw_idx_S2,:); %
          S2_resp_prob_bW= S2_resp_prob(bw_idx_S2,:); %
          touch_types_S2_bW= touch_types_within_S2(bw_idx_S2-max(s1i),:,a);
          
          mean_each_type_BW (a,:,1)= dot(S1_resp_mean_bW, touch_types_S1_bW)/length(touch_types_S1_bW); %S1
          mean_each_type_BW(a,:,2)= dot(S2_resp_mean_bW, touch_types_S2_bW)/length(touch_types_S2_bW); %S2
          
          prob_each_type_BW (a,:,1)= dot(S1_resp_prob_bW, touch_types_S1_bW)/length(touch_types_S1_bW); %S1
          prob_each_type_BW(a,:,2)= dot(S2_resp_prob_bW, touch_types_S2_bW)/length(touch_types_S2_bW); %S2
          
        end 
        
        %USW_mean (a,:,1)=  %USW mean S1, 1 row for each animal
        %USW_mean (a,:,2)= nanmean(mean_S2_to_type (:,a, 1)); % USW mean S2
        
        means_to_plot_S1(a,:,1)= nanmean(mean_S1_to_type (:,a, 1)); %USW, S1
        means_to_plot_S1(a,:,2)= mean(mean_each_type_BW(a,:,1),2); %bw, S1
        means_to_plot_S1(a,:,3)= mean(mean_each_type_MW(a,:,1),2); %mw, S1
        means_to_plot_S2(a,:,1)= nanmean(mean_S2_to_type (:,a, 1)); % USW mean S2
        means_to_plot_S2(a,:,2)= mean(mean_each_type_BW(a,:,2),2); %bw, S2
        means_to_plot_S2(a,:,3)= mean(mean_each_type_MW(a,:,2),2); %mw, S2
        
        probs_to_plot_S1(a,:,1)= nanmean(prob_S1_to_type (:,a, 1)); %USW, S1
        probs_to_plot_S1(a,:,2)= mean(prob_each_type_BW(a,:,1),2); %bw, S1
        probs_to_plot_S1(a,:,3)= mean(prob_each_type_MW(a,:,1),2); %mw, S1
        probs_to_plot_S2(a,:,1)= nanmean(prob_S2_to_type (:,a, 1)); % USW mean S2
        probs_to_plot_S2(a,:,2)= mean(prob_each_type_BW(a,:,2),2); %bw, S2
        probs_to_plot_S2(a,:,3)= mean(prob_each_type_MW(a,:,2),2); %mw, S2

    end
    
    
    
    cell_type = {'uSW','bSW','MW'};
    for cti=1:3 % cell type loop
        %v_s1 = squeeze(frac_all_cells_by_type_and_area(:,1,cti));
        %v_s2 = squeeze(frac_all_cells_by_type_and_area(:,2,cti));
        %[h pv] = ttest(v_s1,v_s2);
        %pv = signrank(v_s1,v_s2);
        %disp(sprintf('overall fraction for %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    %cell_type{cti}, nanmean(v_s1), nanstd(v_s1), nanmean(v_s2), nanstd(v_s2), length(find(~isnan(v_s1))), pv)); 
                
        v_s1 = squeeze(frac_touch_cells_by_type_and_area(:,1,cti));
        v_s2 = squeeze(frac_touch_cells_by_type_and_area(:,2,cti));
        
        
        v_s1_all = squeeze(frac_all_cells_by_type_and_area(:,1,cti));
        v_s2_all = squeeze(frac_all_cells_by_type_and_area(:,2,cti));
        
        [h pv] = ttest(v_s1,v_s2);
        pv = signrank(v_s1,v_s2);
        disp(sprintf('touch only fraction for superfic. %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    cell_type{cti}, nanmean(v_s1), nanstd(v_s1), nanmean(v_s2), nanstd(v_s2), length(find(~isnan(v_s1))), pv));
        
        [h pv] = ttest(v_s1_all,v_s2_all);
        pv = signrank(v_s1_all,v_s2_all);
        disp(sprintf('all cell fraction for superfic. %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    cell_type{cti}, nanmean(v_s1_all), nanstd(v_s1_all), nanmean(v_s2_all), nanstd(v_s2_all), length(find(~isnan(v_s1_all))), pv));
        % plotting barplot of fraction of cells in S1 and S2 with specific
        % type
        x1= [1.5 3.5 5.5]; x2= [2.5 4.5 6.5];   
        color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
        bar(frac_by_type_and_area_ax,x1(cti), nanmean(v_s1),'FaceColor',color_arr(cti,:)); hold on;
        plot (frac_by_type_and_area_ax,x1(cti), v_s1(~isnan(v_s1)),'ko')
        bar(frac_by_type_and_area_ax,x2(cti), nanmean(v_s2),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4);
        xlim= ([0 7]) ; frac_by_type_and_area_ax.YLim= ([0 0.9]); title (frac_by_type_and_area_ax, 'Fraction of touch cells of each type');
        frac_by_type_and_area_ax.TickDir='out';frac_by_type_and_area_ax.Box='off';
        plot (frac_by_type_and_area_ax,x2(cti), v_s2(~isnan(v_s2)),'ko')  
        
        bar (frac_all_ax, x1(cti), nanmean(v_s1_all),'FaceColor',color_arr(cti,:)); hold on;
        plot (frac_all_ax,x1(cti), v_s1_all(~isnan(v_s1_all)),'ko')
        bar(frac_all_ax,x2(cti), nanmean(v_s2_all),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4);
        xlim= ([0 7]) ; frac_all_ax.YLim= ([0 0.08]); title (frac_all_ax, 'Fraction of ALL cells of each type');
        frac_all_ax.TickDir='out';frac_all_ax.Box='off';
        plot (frac_all_ax,x2(cti), v_s2_all(~isnan(v_s2_all)),'ko')  
       
        
%         % plot mean response across types for each 
%         
       % plotting barpolotof mean resposne
%         %type groups in s1 and s2 
          title (ax_mean, 'mean response for cells of each type: averaged across all types')
          bar (ax_mean, x1(cti),mean(nonzeros(means_to_plot_S1(:,:,cti))),'FaceColor',color_arr(cti,:)); hold on; %s1
          plot (ax_mean, x1(cti),(nonzeros(means_to_plot_S1(:,:,cti))),'ko')
          bar (ax_mean, x2(cti),mean(nonzeros(means_to_plot_S2(:,:,cti))),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
          plot (ax_mean, x2(cti),(nonzeros(means_to_plot_S2(:,:,cti))),'ko')
          ax_mean.TickDir='out';ax_mean.Box='off';
          
%         % plotting barpolot of mean prob of respon;
%         %type groups in s1 and s2 
          title (ax_2, 'prob response for cells of each type: averaged across all types')
          bar (ax_2, x1(cti),mean(nonzeros(probs_to_plot_S1(:,:,cti))),'FaceColor',color_arr(cti,:)); hold on; %s1
          plot (ax_2, x1(cti),(nonzeros(probs_to_plot_S1(:,:,cti))),'ko')
          bar (ax_2, x2(cti),mean(nonzeros(probs_to_plot_S2(:,:,cti))),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
          plot (ax_2, x2(cti),(nonzeros(probs_to_plot_S2(:,:,cti))),'ko')
          ax_2.TickDir='out';ax_2.Box='off'; 
       
          A = nonzeros(means_to_plot_S1(:,:,cti));
          B = nonzeros(means_to_plot_S2(:,:,cti));
          [h,pv]=ttest(A,B);
          pval(cti)= pv;
      
%       
    end 
   
end 

function figure_compare_s1_s2_touch_deep
   all_dat = get_s1s2_all_dat;
    anims_used = {'an014351','an016623','an017517','an016650','an016652','an014363','an017510'};
    ani_used = find(ismember(all_dat.anims, anims_used));
    
    %find(all_dat.has_deep_data_s1s2); 
    % --- make the blank plots
    fh = figure('Position', [0 0 2000 500]); 
    
    frac_by_type_and_area_ax = subplot('Position', [.05 .05 .15 .8]);  
    hold on;
    ax_2 = subplot('Position', [.25 .05 .15 .8]);  
    hold on;
    ax_mean= subplot('Position', [.45 .05 .15 .8]);   
    hold on;
    frac_all_ax = subplot('Position', [.65 .05 .15 .8]);  
    hold on;
   
    
    % --- fraction of cells by area and type
    frac_all_cells_by_type_and_area = nan*zeros(length(all_dat.anims), 2, 3); % animal X s1/s2 X uSW/bSW/MW
    frac_touch_cells_by_type_and_area = nan*zeros(length(all_dat.anims), 2, 3); % animal X s1/s2 X uSW/bSW/MW
    mean_S1_to_type= nan*zeros(4,length(all_dat.anims),3);
    mean_S2_to_type= nan*zeros(4,length(all_dat.anims),3);
    prob_S1_to_type= nan*zeros(4,length(all_dat.anims),3);
    prob_S2_to_type= nan*zeros(4,length(all_dat.anims),3);
    %mean_prob_resp_w1p= nan*zeros(length(all_dat.anims),2,3); %animal, s1/s2, usw/bsw/mw -- this will contain mean probability of response to w1p for each cell type in s1 and s2
    %S1_MW_resp =nan*zeros(100,4,length(all_dat.anims))
    
    d_s1 = 2; % most superficial of merged subvolumes
    d_s2 = 4; 
  
    
    for ai=1:length(ani_used)
        a = ani_used(ai);

        % figure out which neurons go where
%        s1i = find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s1_id_range{a}(2));
%        s2i = find(all_dat.anim_data{a}.ids > all_dat.s2_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s2_id_range{a}(2));
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);
        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);

        % which date should we use? second non-compound date
        non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);

        d_s1_used = d_s1;
        d_s2_used = d_s2; 
        
        
        % animal-specific stuff
        if (strcmp(all_dat.anims{a}, 'an014359')) %has no merging, so need to use the best day, but also need to seperate superficial and deep!
            d_s2_used = 1;
            %s1i= s1i(all_dat.anim_data{3}.ids(s1i)<20000000);
%            s1i= s1i(all_dat.anim_data{3}.ids(s1i)>20000000);
            %s2i= s2i(all_dat.anim_data{3}.ids(s2i)<920000000);
%            s2i= s2i(all_dat.anim_data{3}.ids(s2i)>920000000);
        elseif (strcmp(all_dat.anims{a}, 'an016652'))
            d_s1_used = 4;
            d_s2_used = 2;
            
        end
        
       
        frac_all_cells_by_type_and_area (a, 1, 1) = length(intersect(s1i,all_dat.types_by_idx{a}.usw_by_day{d_s1_used}))/length(s1i);
        frac_all_cells_by_type_and_area (a, 2, 1) = length(intersect(s2i,all_dat.types_by_idx{a}.usw_by_day{d_s2_used}))/length(s2i);

        frac_all_cells_by_type_and_area (a, 1, 2) = length(intersect(s1i,all_dat.types_by_idx{a}.bsw_by_day{d_s1_used}))/length(s1i);
        frac_all_cells_by_type_and_area (a, 2, 2) = length(intersect(s2i,all_dat.types_by_idx{a}.bsw_by_day{d_s2_used}))/length(s2i);

        frac_all_cells_by_type_and_area (a, 1, 3) = length(intersect(s1i,all_dat.types_by_idx{a}.mw_by_day{d_s1_used}))/length(s1i);
        frac_all_cells_by_type_and_area (a, 2, 3) = length(intersect(s2i,all_dat.types_by_idx{a}.mw_by_day{d_s2_used}))/length(s2i);

        all_touch_s1_idx = [all_dat.types_by_idx{a}.usw_by_day{d_s1_used} all_dat.types_by_idx{a}.bsw_by_day{d_s1_used} all_dat.types_by_idx{a}.mw_by_day{d_s1_used}];
        all_touch_s2_idx = [all_dat.types_by_idx{a}.usw_by_day{d_s2_used} all_dat.types_by_idx{a}.bsw_by_day{d_s2_used} all_dat.types_by_idx{a}.mw_by_day{d_s2_used}];
       
        frac_touch_cells_by_type_and_area (a, 1, 1) = length(intersect(s1i,all_dat.types_by_idx{a}.usw_by_day{d_s1_used}))/length(all_touch_s1_idx);
        frac_touch_cells_by_type_and_area (a, 2, 1) = length(intersect(s2i,all_dat.types_by_idx{a}.usw_by_day{d_s2_used}))/length(all_touch_s2_idx);

        frac_touch_cells_by_type_and_area (a, 1, 2) = length(intersect(s1i,all_dat.types_by_idx{a}.bsw_by_day{d_s1_used}))/length(all_touch_s1_idx);
        frac_touch_cells_by_type_and_area (a, 2, 2) = length(intersect(s2i,all_dat.types_by_idx{a}.bsw_by_day{d_s2_used}))/length(all_touch_s2_idx);

        frac_touch_cells_by_type_and_area (a, 1, 3) = length(intersect(s1i,all_dat.types_by_idx{a}.mw_by_day{d_s1_used}))/length(all_touch_s1_idx);
        frac_touch_cells_by_type_and_area (a, 2, 3) = length(intersect(s2i,all_dat.types_by_idx{a}.mw_by_day{d_s2_used}))/length(all_touch_s2_idx);
            
        % make matrix of MW Responsive cells
            mw_idx_S1= intersect(all_dat.types_by_idx{a}.mw_by_day{d_s1_used},s1i); %find ids of multiwhisker cells in each subvolume
            mw_idx_S2= intersect(all_dat.types_by_idx{a}.mw_by_day{d_s2_used},s2i);
            bw_idx_S1= intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s1_used},s1i);
            bw_idx_S2= intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s2_used},s2i);

            w1p_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d_s1_used}, s1i); %find indexes of each type of touch in s1 and s2
            w1r_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d_s1_used}, s1i);
            w2p_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d_s1_used}, s1i);
            w2r_idx_S1= intersect(all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d_s1_used}, s1i);    

            w1p_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d_s2_used}, s2i);
            w1r_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d_s2_used}, s2i);
            w2p_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d_s2_used}, s2i);
            w2r_idx_S2= intersect(all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d_s2_used}, s2i);  

            for i= 1: length(s1i) %loop thorugh the possible all cell ids and make a matrix of cell id x type (w1p, w1r, w2p, w2r) per animal in S1 
               touch_types_within_S1(i,1, a)= any(w1p_idx_S1(:)== s1i(i));
               touch_types_within_S1(i,2,a)= any(w1r_idx_S1(:)== s1i(i));
               touch_types_within_S1(i,3,a)= any(w2p_idx_S1(:)== s1i(i));
               touch_types_within_S1(i, 4,a)= any(w2r_idx_S1(:)== s1i(i));   
            end      

             for i= 1: length(s2i) % same for S2
               touch_types_within_S2(i,1, a)= any(w1p_idx_S2(:)== s2i(i));
               touch_types_within_S2(i,2,a)= any(w1r_idx_S2(:)== s2i(i));
               touch_types_within_S2(i,3,a)= any(w2p_idx_S2(:)== s2i(i));
               touch_types_within_S2(i, 4,a)= any(w2r_idx_S2(:)== s2i(i));   
             end 
            
        % so touch types within S1 gives for each cell whether it was
        % response to w1p, then w1r, then w2p, then w2r-- full matrix, can
        % subselect by subtype of touch later. 
        
        S1_resp_mean= nan*zeros(max(s1i),4);
        S2_resp_mean= nan*zeros(max(s2i),4);
        S1_resp_prob= nan*zeros(max(s1i),4);
        S2_resp_prob= nan*zeros(max(s2i),4);
        %responses to each touch type (mean)
        for tt= 1:4 % loop through 4 touch types: w1p, w1r, w2p, w2r
           if tt==1 response=all_dat.types_by_idx{a}.ev_touch_w1p_by_day, mean_use=all_dat.anim_mats{a}.meanRespW1p; prob_use=all_dat.anim_mats{a}.probRespW1p; end 
           if tt==3 response=all_dat.types_by_idx{a}.ev_touch_w2p_by_day; mean_use= all_dat.anim_mats{a}.meanRespW2p; prob_use=all_dat.anim_mats{a}.probRespW2p;end 
           if tt==2 response=all_dat.types_by_idx{a}.ev_touch_w1r_by_day; mean_use= all_dat.anim_mats{a}.meanRespW1r;prob_use=all_dat.anim_mats{a}.probRespW1r;end 
           if tt==4 response=all_dat.types_by_idx{a}.ev_touch_w2r_by_day; mean_use= all_dat.anim_mats{a}.meanRespW2r;prob_use=all_dat.anim_mats{a}.probRespW2r;end 
          
           %example: mean_S1_to_type(1,:,1) is
           %the mean response in S1 for all animals to of usw cells to w1p.
           %this works for USw but not hte others, need to do those
           %seperately. 
           %order: usw,bw, mw
           mean_S1_to_type(tt,a,1)= mean(mean_use(d_s1_used, intersect(response{d_s1_used},all_dat.types_by_idx{a}.usw_by_day{d_s1_used})));
           
           mean_S2_to_type(tt, a, 1)= mean(mean_use(d_s2_used, intersect(response{d_s2_used},all_dat.types_by_idx{a}.usw_by_day{d_s2_used})));
           
           prob_S1_to_type(tt,a,1)= mean(prob_use(d_s1_used, intersect(response{d_s1_used},all_dat.types_by_idx{a}.usw_by_day{d_s1_used})));
           
           prob_S2_to_type(tt, a, 1)= mean(prob_use(d_s2_used, intersect(response{d_s2_used},all_dat.types_by_idx{a}.usw_by_day{d_s2_used})));
           
          % for each animal, construct a matrix of mean response of each
          % neuron to each touch type. we can take the dot product of this
          % matrix with the matrix listing what is resposive to what 
         
          S1_resp_mean (min(s1i):max(s1i),tt)= mean_use(d_s1_used, s1i);
          S2_resp_mean (min(s2i):max(s2i),tt)= mean_use(d_s2_used, s2i);
          S1_resp_prob (min(s1i):max(s1i),tt)= prob_use(d_s1_used, s1i);
          S2_resp_prob (min(s2i):max(s2i),tt)= prob_use(d_s2_used, s2i);
          
          S1_resp_mean_MW= S1_resp_mean(mw_idx_S1,:); %use only MW cells for dot prod.
          S1_resp_prob_MW= S1_resp_prob(mw_idx_S1,:); %use only MW cells for dot prod.

          touch_types_S1_MW= touch_types_within_S1(mw_idx_S1, :,a); %use only MW cells for dot prod
          
          S2_resp_mean_MW= S2_resp_mean(mw_idx_S2,:); %
          S2_resp_prob_MW= S2_resp_prob(mw_idx_S2,:); %

          touch_types_S2_MW= touch_types_within_S2(mw_idx_S2-max(s1i),:,a);
          
         
          mean_each_type_MW (a,:,1)= dot(S1_resp_mean_MW, touch_types_S1_MW)/length(touch_types_S1_MW); %S1
          mean_each_type_MW(a,:,2)= dot(S2_resp_mean_MW, touch_types_S2_MW)/length(touch_types_S2_MW); %S2
          
          prob_each_type_MW (a,:,1)= dot(S1_resp_prob_MW, touch_types_S1_MW)/length(touch_types_S1_MW); %S1
          prob_each_type_MW(a,:,2)= dot(S2_resp_prob_MW, touch_types_S2_MW)/length(touch_types_S2_MW); %S2
          
          %ok now do same dot product thing for bidir cells
          S1_resp_mean_bW= S1_resp_mean(bw_idx_S1,:); %use only MW cells for dot prod.
          S1_resp_prob_bW= S1_resp_prob(bw_idx_S1,:);
          
          touch_types_S1_bW= touch_types_within_S1(bw_idx_S1, :,a); %use only MW cells for dot prod
          S2_resp_mean_bW= S2_resp_mean(bw_idx_S2,:); %
          S2_resp_prob_bW= S2_resp_prob(bw_idx_S2,:); %
          touch_types_S2_bW= touch_types_within_S2(bw_idx_S2-max(s1i),:,a);
          
          mean_each_type_BW (a,:,1)= dot(S1_resp_mean_bW, touch_types_S1_bW)/length(touch_types_S1_bW); %S1
          mean_each_type_BW(a,:,2)= dot(S2_resp_mean_bW, touch_types_S2_bW)/length(touch_types_S2_bW); %S2
          
          prob_each_type_BW (a,:,1)= dot(S1_resp_prob_bW, touch_types_S1_bW)/length(touch_types_S1_bW); %S1
          prob_each_type_BW(a,:,2)= dot(S2_resp_prob_bW, touch_types_S2_bW)/length(touch_types_S2_bW); %S2
          
        end 
        
        %USW_mean (a,:,1)=  %USW mean S1, 1 row for each animal
        %USW_mean (a,:,2)= nanmean(mean_S2_to_type (:,a, 1)); % USW mean S2
        
        means_to_plot_S1(a,:,1)= nanmean(mean_S1_to_type (:,a, 1)); %USW, S1
        means_to_plot_S1(a,:,2)= mean(mean_each_type_BW(a,:,1),2); %bw, S1
        means_to_plot_S1(a,:,3)= mean(mean_each_type_MW(a,:,1),2); %mw, S1
        means_to_plot_S2(a,:,1)= nanmean(mean_S2_to_type (:,a, 1)); % USW mean S2
        means_to_plot_S2(a,:,2)= mean(mean_each_type_BW(a,:,2),2); %bw, S2
        means_to_plot_S2(a,:,3)= mean(mean_each_type_MW(a,:,2),2); %mw, S2
        
        probs_to_plot_S1(a,:,1)= nanmean(prob_S1_to_type (:,a, 1)); %USW, S1
        probs_to_plot_S1(a,:,2)= mean(prob_each_type_BW(a,:,1),2); %bw, S1
        probs_to_plot_S1(a,:,3)= mean(prob_each_type_MW(a,:,1),2); %mw, S1
        probs_to_plot_S2(a,:,1)= nanmean(prob_S2_to_type (:,a, 1)); % USW mean S2
        probs_to_plot_S2(a,:,2)= mean(prob_each_type_BW(a,:,2),2); %bw, S2
        probs_to_plot_S2(a,:,3)= mean(prob_each_type_MW(a,:,2),2); %mw, S2

    end
    
    
    
    cell_type = {'uSW','bSW','MW'};
    for cti=1:3 % cell type loop
        %v_s1 = squeeze(frac_all_cells_by_type_and_area(:,1,cti));
        %v_s2 = squeeze(frac_all_cells_by_type_and_area(:,2,cti));
        %[h pv] = ttest(v_s1,v_s2);
        %pv = signrank(v_s1,v_s2);
        %disp(sprintf('overall fraction for %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    %cell_type{cti}, nanmean(v_s1), nanstd(v_s1), nanmean(v_s2), nanstd(v_s2), length(find(~isnan(v_s1))), pv)); 
                
        v_s1 = squeeze(frac_touch_cells_by_type_and_area(:,1,cti));
        v_s2 = squeeze(frac_touch_cells_by_type_and_area(:,2,cti));
        
        
        v_s1_all = squeeze(frac_all_cells_by_type_and_area(:,1,cti));
        v_s2_all = squeeze(frac_all_cells_by_type_and_area(:,2,cti));
        
        [h pv] = ttest(v_s1,v_s2);
        pv = signrank(v_s1,v_s2);
        disp(sprintf('touch only fraction for superfic. %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    cell_type{cti}, nanmean(v_s1), nanstd(v_s1), nanmean(v_s2), nanstd(v_s2), length(find(~isnan(v_s1))), pv));
        
        [h pv] = ttest(v_s1_all,v_s2_all);
        pv = signrank(v_s1_all,v_s2_all);
        disp(sprintf('all cell fraction for superfic. %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    cell_type{cti}, nanmean(v_s1_all), nanstd(v_s1_all), nanmean(v_s2_all), nanstd(v_s2_all), length(find(~isnan(v_s1_all))), pv));
        % plotting barplot of fraction of cells in S1 and S2 with specific
        % type
        x1= [1.5 3.5 5.5]; x2= [2.5 4.5 6.5];   
        color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
        bar(frac_by_type_and_area_ax,x1(cti), nanmean(v_s1),'FaceColor',color_arr(cti,:)); hold on;
        plot (frac_by_type_and_area_ax,x1(cti), v_s1(~isnan(v_s1)),'ko')
        bar(frac_by_type_and_area_ax,x2(cti), nanmean(v_s2),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4);
        xlim= ([0 7]) ; frac_by_type_and_area_ax.YLim= ([0 0.9]); title (frac_by_type_and_area_ax, 'Fraction of touch cells of each type');
        frac_by_type_and_area_ax.TickDir='out';frac_by_type_and_area_ax.Box='off';
        plot (frac_by_type_and_area_ax,x2(cti), v_s2(~isnan(v_s2)),'ko')  
        
        bar (frac_all_ax, x1(cti), nanmean(v_s1_all),'FaceColor',color_arr(cti,:)); hold on;
        plot (frac_all_ax,x1(cti), v_s1_all(~isnan(v_s1_all)),'ko')
        bar(frac_all_ax,x2(cti), nanmean(v_s2_all),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4);
        xlim= ([0 7]) ; frac_all_ax.YLim= ([0 0.08]); title (frac_all_ax, 'Fraction of ALL cells of each type');
        frac_all_ax.TickDir='out';frac_all_ax.Box='off';
        plot (frac_all_ax,x2(cti), v_s2_all(~isnan(v_s2_all)),'ko')  
       
        
%         % plot mean response across types for each 
%         
       % plotting barpolotof mean resposne
%         %type groups in s1 and s2 
          title (ax_mean, 'mean response for cells of each type: averaged across all types')
          bar (ax_mean, x1(cti),nanmean(nonzeros(means_to_plot_S1(:,:,cti))),'FaceColor',color_arr(cti,:)); hold on; %s1
          plot (ax_mean, x1(cti),(nonzeros(means_to_plot_S1(:,:,cti))),'ko')
          bar (ax_mean, x2(cti),nanmean(nonzeros(means_to_plot_S2(:,:,cti))),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
          plot (ax_mean, x2(cti),(nonzeros(means_to_plot_S2(:,:,cti))),'ko')
          ax_mean.TickDir='out';ax_mean.Box='off';
          
%         % plotting barpolot of mean prob of respon;
%         %type groups in s1 and s2 
          title (ax_2, 'prob response for cells of each type: averaged across all types')
          bar (ax_2, x1(cti),nanmean(nonzeros(probs_to_plot_S1(:,:,cti))),'FaceColor',color_arr(cti,:)); hold on; %s1
          plot (ax_2, x1(cti),(nonzeros(probs_to_plot_S1(:,:,cti))),'ko')
          bar (ax_2, x2(cti),nanmean(nonzeros(probs_to_plot_S2(:,:,cti))),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
          plot (ax_2, x2(cti),(nonzeros(probs_to_plot_S2(:,:,cti))),'ko')
          ax_2.TickDir='out';ax_2.Box='off'; 
       
          A = nonzeros(means_to_plot_S1(:,:,cti));
          B = nonzeros(means_to_plot_S2(:,:,cti));
          %[h,pv]=ttest(A,B);
          pval(cti)= pv;
      
%       
    end 
    
end 

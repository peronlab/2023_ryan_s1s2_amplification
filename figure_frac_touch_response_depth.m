   % computes the fraction of touch response in each area per depth. 
   
function frac_touch_response_depth
    all_dat = get_s1s2_all_dat;
    vali= find(all_dat.has_deep_data_s1s2); 
    cts = {'usw', 'bsw', 'mw'}; 
    %ratio_mat = zeros(3,length(vali));
    bins=85:50:435;
    mw_id= zeros(300,2,2); bsw_id= zeros(300,2,2); usw_id= zeros(300,2,2);
    
for v=1:length(vali)
    ai=vali(v);
    
    for ii=1:2 %loop thorugh S1 or S2
        if ii==1 area_nrni = get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, 's1_manual_restrict', all_dat, ai);
            d_s=1;
            d_d=2; 
            if (strcmp(all_dat.anims{ai}, 'an016652'))
               d_s=3;
               d_d=4; 
            end 
        end 
        if ii==2 area_nrni = get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, 's2_manual_restrict', all_dat, ai);
            d_s=3;
            d_d=4;
            if (strcmp(all_dat.anims{ai}, 'an016652')) 
                d_s=1;
                d_d=2;
            end 
            if (strcmp(all_dat.anims{ai}, 'an014359'))
                d_s= 1;
                d_d= 1;
            end 
        end 
    
        for jj=1:2
            if jj==1 dat_dayi= d_s; end %superficial 
            if jj==2 dat_dayi= d_d; end %deep
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
            
            mw_idx= intersect(area_nrni,all_dat.types_by_idx{ai}.mw_by_day{dat_dayi}); 
            bsw_idx= intersect(area_nrni, all_dat.types_by_idx{ai}.bsw_by_day{dat_dayi});
            usw_idx= intersect(area_nrni, all_dat.types_by_idx{ai}.usw_by_day{dat_dayi});

            if ii==1 && jj==1 S1s= mu_resp_per_cell; mw_S1s= mw_idx; bsw_S1s= bsw_idx; usw_S1s= usw_idx; end %we are finding mean response in 4 different subvolumes, by touch type cell actually responds to.
            if ii==1 && jj==2 S1d= mu_resp_per_cell; mw_S1d= mw_idx; bsw_S1d= bsw_idx; usw_S1d= usw_idx; end 
            if ii==2 && jj==1 S2s= mu_resp_per_cell; mw_S2s= mw_idx; bsw_S2s= bsw_idx; usw_S2s= usw_idx; end 
            if ii==2 && jj==2 S2d= mu_resp_per_cell; mw_S2d= mw_idx; bsw_S2d= bsw_idx; usw_S2d= usw_idx; end 
            
            
            %mw_id  (1:length(mw_idx),jj,ii)= mw_idx; % structure is mw cells/super or deep/area S1 or S2 
            %bsw_id (1:length(bsw_idx),jj,ii)= bsw_idx;
            %usw_id (1:length(usw_idx),jj,ii)= usw_idx; 
            
        end 
    end
    for c=1:3
           if c==1 ctid_S1= [usw_S1s' ,usw_S1d']; ctid_S2= [usw_S2s' ,usw_S2d'];end 
           if c==2 ctid_S1= [bsw_S1s' ,bsw_S1d']; ctid_S2= [bsw_S2s' ,bsw_S2d'];  end 
           if c==3 ctid_S1= [mw_S1s' ,mw_S1d']; ctid_S2= [mw_S2s' ,mw_S2d']; end 
           
           %cells_of_type_S1= unique(nonzeros([ctid(:,1,1); ctid(:,2,1)])); % concatenate cells of both subvolumes, still 1 area
           depth_cells_S1=  all_dat.anim_data{ai}.z_um(ctid_S1); % depth of those cells 
           %cells_of_type_S2=unique(nonzeros([ctid(:,1,2); ctid(:,2,2)])); % concatenate cells of both subvolumes, still 1 area
           depth_cells_S2= all_dat.anim_data{ai}.z_um(ctid_S2); % depth of those cells 
           
           for i= 1: length(bins)-1 %set up depth bins
              cells_in_bin_of_type_s1= find((depth_cells_S1>= bins(i))&(depth_cells_S1 <= bins(i+1)));
              idx_cells_in_bin_type_s1= ctid_S1(cells_in_bin_of_type_s1);
              cells_in_bin_of_type_s2= find((depth_cells_S2>= bins(i))&(depth_cells_S2 <= bins(i+1)));
              idx_cells_in_bin_type_s2= ctid_S2(cells_in_bin_of_type_s2);
              %accoutn for both superficial and deep by calculating means
              %for both and if there are no cells in this bin in that
              %catergory, itll just not add to the sum
              means_s1(c,i)= nansum([nansum(S1s(idx_cells_in_bin_type_s1)) nansum(S1d(idx_cells_in_bin_type_s1))]);
              
              means_s2(c,i)= nansum([nansum(S2s(idx_cells_in_bin_type_s2)) nansum(S2d(idx_cells_in_bin_type_s2))]);
      
           end 
    
     frac_resp_bin (c,:,v)= means_s1(c,:)/nansum(means_s1(c,:));   
     frac_resp_bin_s2 (c,:,v)= means_s2(c,:)/nansum(means_s2(c,:));
    
    end
    sum_an_S1= sum(means_s1);
    sum_an_S2= sum(means_s2); 
    mean_S1_all (:,:,v)= means_s1./sum_an_S1;
    mean_S2_all (:,:,v)= means_s2./sum_an_S2;
end 

plot_frac_1= mean(frac_resp_bin,3);
plot_frac_2= mean(frac_resp_bin_s2,3); 
%mean_S2_all= (mean_S2_all)
plot_frac_all_1= nanmean(mean_S1_all,3);
plot_frac_all_2= nanmean(mean_S2_all,3);

    %%  %%
    fh = figure('Position', [0 0 1000 500]); 
    frac_at_d_S1 = subplot('Position', [.05 .1 0.2 .8]);  
    hold(frac_at_d_S1,'on');
    frac_at_d_S2 = subplot('Position', [.3 .1 .2 .8]);  
    hold(frac_at_d_S2,'on');
    
        color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
        bins_to_plot= (bins(:,1:end-1) + bins(:,2:end)) / 2; %Find midpoint of bin to plot in 
        
        for yy= 1:3
            plot (frac_at_d_S1,plot_frac_all_1(yy,:), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)
            plot (frac_at_d_S2,plot_frac_all_2(yy,:), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)

            %plot (norm_frac_at_d_S1,depth_fracs_S1_overall(yy,:,1)/max(depth_fracs_S1_overall(yy,:,1)), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)
            %plot (norm_frac_at_d_S2,depth_fracs_S2_overall(yy,:,1)/max(depth_fracs_S2_overall(yy,:,1)), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)
        end    

        set(frac_at_d_S1, 'YDir','reverse')
        frac_at_d_S1.TickDir='out';
        %ylim ([100 400]); hold on
        xlabel (frac_at_d_S1, 'Fraction of touch response per type')
        ylabel (frac_at_d_S1, 'Depth')
        title (frac_at_d_S1,'s1')
        
        
        set(frac_at_d_S2, 'YDir','reverse')
        frac_at_d_S2.TickDir='out';
        xlabel (frac_at_d_S2, 'Fraction of touch response per type');
        ylabel (frac_at_d_S2, 'Depth')
        title (frac_at_d_S2,'s2');
        %ylim( [100 400])
        %aa = axis(norm_frac_at_d_S1) ; axis(norm_frac_at_d_S1, [-0.1 1.1 aa(3) aa(4)]);

end
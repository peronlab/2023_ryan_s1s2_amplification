% depth plots 
% 1: find all touch cell ids in S1 
% 2:find all MW touch cell ids in S1 
% 3: make depth bins, and count #mw/#touch
% 1: find multiwhisker ids for each animal 
% 2: determine what time bins those go in 
% 3: plot 

% for this, incorporating both subvolumes by just concatenating arrays for
% the 2 subvols of each area.
%% 
%%
function figure_depth
    all_dat = get_s1s2_all_dat;
    ani_used = find(all_dat.has_deep_data_s1s2);
    
%    bins= [50 100 150 200 250 300 350 400 450]; 
%    bins = 50:50:450;
    bins= 75:50:425;

    frac_at_depth_S1= nan*zeros(3, length(bins)-1, length(ani_used)); %  uSW/bSW/MW x depth bin x animal
    frac_at_depth_S2= nan*zeros(3, length(bins)-1, length(ani_used)); %  uSW/bSW/MW x depth bin x animal

    % --- make the blank plots
    fh = figure('Position', [0 0 1000 500]); 
    frac_at_d_S1 = subplot('Position', [.05 .1 0.2 .8]);  
    hold(frac_at_d_S1,'on');
    frac_at_d_S2 = subplot('Position', [.3 .1 .2 .8]);  
    hold(frac_at_d_S2,'on');
    
    norm_frac_at_d_S1 = subplot('Position', [.55 .1 0.2 .8]);  
    hold(norm_frac_at_d_S1,'on');
    norm_frac_at_d_S2 = subplot('Position', [.8 .1 .2 .8]);  
    hold(norm_frac_at_d_S2,'on');

    d_s1_s = 1; % most superficial of merged subvolumes
    d_s1_d= 2;
    d_s2_s = 3; 
    d_s2_d= 4; 
    %% 
    for ai= 1:length(ani_used)
        a = ani_used(ai);

        % figure out which neurons go where
        %s1i = find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s1_id_range{a}(2));
        %s2i = find(all_dat.anim_data{a}.ids > all_dat.s2_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s2_id_range{a}(2));

        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);

        % which date should we use? second non-compound date
        non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);

        d_s1_used_s= d_s1_s;
        d_s1_used_d= d_s1_d;
        d_s2_used_s= d_s2_s;
        d_s2_used_d= d_s2_d;
        
        % animal-specific stuff
        if (strcmp(all_dat.anims{a}, 'an014359'))
            d_s2_used_s = 1;
            d_s2_used_d= 1;
        elseif (strcmp(all_dat.anims{a}, 'an016652'))
            d_s1_used_s = 3;
            d_s1_used_d= 4;
            d_s2_used_s = 1;
            d_s2_used_d= 2;
        end
%%
        all_touch_s1_idx = s1i;%[all_dat.types_by_idx{a}.usw_by_day{d_s1_used_s} all_dat.types_by_idx{a}.usw_by_day{d_s1_used_d} all_dat.types_by_idx{a}.bsw_by_day{d_s1_used_s} all_dat.types_by_idx{a}.bsw_by_day{d_s1_used_d}  all_dat.types_by_idx{a}.mw_by_day{d_s1_used_s} all_dat.types_by_idx{a}.mw_by_day{d_s1_used_d}]; 
        mw_idx_S1= [intersect(all_dat.types_by_idx{a}.mw_by_day{d_s1_used_s},s1i); intersect(all_dat.types_by_idx{a}.mw_by_day{d_s1_used_d},s1i)]; %find ids of multiwhisker cells in each subvolume
        bsw_idx_S1= [intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s1_used_s},s1i);intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s1_used_d},s1i)];
        usw_idx_S1= [intersect(all_dat.types_by_idx{a}.usw_by_day{d_s1_used_s},s1i);intersect(all_dat.types_by_idx{a}.usw_by_day{d_s1_used_d},s1i)];

        % S2 
        all_touch_s2_idx = s2i;%[all_dat.types_by_idx{a}.usw_by_day{d_s2_used_s} all_dat.types_by_idx{a}.usw_by_day{d_s2_used_d} all_dat.types_by_idx{a}.bsw_by_day{d_s2_used_s} all_dat.types_by_idx{a}.bsw_by_day{d_s2_used_d}  all_dat.types_by_idx{a}.mw_by_day{d_s2_used_s} all_dat.types_by_idx{a}.mw_by_day{d_s2_used_d}]; 
        mw_idx_S2= [intersect(all_dat.types_by_idx{a}.mw_by_day{d_s2_used_s},s2i); intersect(all_dat.types_by_idx{a}.mw_by_day{d_s2_used_d},s2i)]; %find ids of multiwhisker cells in each subvolume
        bsw_idx_S2= [intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s2_used_s},s2i);intersect(all_dat.types_by_idx{a}.bsw_by_day{d_s2_used_d},s2i)];
        usw_idx_S2= [intersect(all_dat.types_by_idx{a}.usw_by_day{d_s2_used_s},s2i);intersect(all_dat.types_by_idx{a}.usw_by_day{d_s2_used_d},s2i)];
    
        depth_all_s1 = all_dat.anim_data{a}.z_um(s1i); %depth of all touch neurons 
        depth_s1= all_dat.anim_data{a}.z_um(all_touch_s1_idx); %depth of all touch neurons 
        depth_mW_s1= all_dat.anim_data{a}.z_um(mw_idx_S1);
        depth_bsW_s1= all_dat.anim_data{a}.z_um(bsw_idx_S1);
        depth_usw_s1= all_dat.anim_data{a}.z_um(usw_idx_S1);
        
        % depths in S2 
        depth_all_s2 = all_dat.anim_data{a}.z_um(s2i); %depth of all touch neurons 
        depth_S2= all_dat.anim_data{a}.z_um(all_touch_s2_idx); %depth of all touch neurons 
        depth_mW_s2= all_dat.anim_data{a}.z_um(mw_idx_S2);
        depth_bsW_s2= all_dat.anim_data{a}.z_um(bsw_idx_S2);
        depth_usw_s2= all_dat.anim_data{a}.z_um(usw_idx_S2);
        
        %bins= [50 100 150 200 250 300 350 400 450]; 
       %bins=linspace (0,400,7);
        for i= 1:length(bins)-1
%            norm_s1 = length(find((depth_s1>= bins(i))&(depth_s1<= bins(i+1)))); % this is touch cells
            norm_s1 = length(find((depth_all_s1>= bins(i))&(depth_all_s1<= bins(i+1)))); % all cells
            frac_at_depth_S1(1,i,ai)= length(find((depth_usw_s1 >= bins(i))&(depth_usw_s1 <= bins(i+1))))/norm_s1;
            frac_at_depth_S1(3,i,ai)= length(find((depth_mW_s1 >= bins(i))&(depth_mW_s1 <= bins(i+1))))/norm_s1;
            frac_at_depth_S1(2,i,ai)= length(find((depth_bsW_s1 >= bins(i))&(depth_bsW_s1 <= bins(i+1))))/norm_s1;
            
%            norm_s2 = length(find((depth_s2>= bins(i))&(depth_s2<= bins(i+1)))); % this is touch cells
            norm_s2 = length(find((depth_all_s2>= bins(i))&(depth_all_s2<= bins(i+1)))); % all cells
            frac_at_depth_S2(1,i,ai)= length(find((depth_usw_s2 >= bins(i))&(depth_usw_s2 <= bins(i+1))))/norm_s2;
            frac_at_depth_S2(3,i,ai)= length(find((depth_mW_s2 >= bins(i))&(depth_mW_s2 <= bins(i+1))))/norm_s2;
            frac_at_depth_S2(2,i,ai)= length(find((depth_bsW_s2 >= bins(i))&(depth_bsW_s2 <= bins(i+1))))/norm_s2;
        end
    end 
    
    for zz=1: 3 %make a matrix of depth fracs computed across animals size usw/bsw/mw x depth bin
        mat_temp=frac_at_depth_S1(zz,:,:);
        mat_temp2= frac_at_depth_S2(zz,:,:);
        depth_fracs_S1_overall(zz,:) = nanmean(mat_temp,3); % take the mean across animals, in each depth bin
        depth_fracs_S2_overall (zz, :)= nanmean(mat_temp2,3);
    end 
  
   %%
        color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
        bins_to_plot= (bins(:,1:end-1) + bins(:,2:end)) / 2; %Find midpoint of bin to plot in 
        
        for yy= 1:3
            plot (frac_at_d_S1,depth_fracs_S1_overall(yy,:,1), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)
            plot (frac_at_d_S2,depth_fracs_S2_overall(yy,:,1), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)

            plot (norm_frac_at_d_S1,depth_fracs_S1_overall(yy,:,1)/max(depth_fracs_S1_overall(yy,:,1)), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)
            plot (norm_frac_at_d_S2,depth_fracs_S2_overall(yy,:,1)/max(depth_fracs_S2_overall(yy,:,1)), bins_to_plot, 'color', color_arr(yy,:),'LineWidth', 6)
        end    

        set(frac_at_d_S1, 'YDir','reverse')
        frac_at_d_S1.TickDir='out';
        xlabel (frac_at_d_S1, 'Fraction of touch neurons with type')
        ylabel (frac_at_d_S1, 'Depth')
        title (frac_at_d_S1,'s1')
 

        set(norm_frac_at_d_S1, 'YDir','reverse')
        norm_frac_at_d_S1.TickDir='out';
        xlabel (norm_frac_at_d_S1, 'Peak-normalized');
        ylabel (norm_frac_at_d_S1, 'Depth')
        title (norm_frac_at_d_S1,'s1');
        aa = axis(norm_frac_at_d_S1) ; axis(norm_frac_at_d_S1, [-0.1 1.1 aa(3) aa(4)]);


        set(frac_at_d_S2, 'YDir','reverse')
        frac_at_d_S2.TickDir='out';
        xlabel (frac_at_d_S2, 'Fraction of touch neurons with type');
        ylabel (frac_at_d_S2, 'Depth')
        title (frac_at_d_S2,'s2')


        set(norm_frac_at_d_S2, 'YDir','reverse')
        norm_frac_at_d_S2.TickDir='out';
        xlabel (norm_frac_at_d_S2, 'Peak-normalized');
        ylabel (norm_frac_at_d_S2, 'Depth')
        title (norm_frac_at_d_S2,'s2')
        aa = axis(norm_frac_at_d_S2) ; axis(norm_frac_at_d_S2, [-0.1 1.1 aa(3) aa(4)]);

end   
    

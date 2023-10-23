%
% Figure for following panels:
%   - fraction of all cells by area that have specific type (uSW, bSW, MW)
%   - mean encoding score for touch cell (may want to take top 5% and mean of that, for example)
%   - ???
%
function figure_compare_s1_s2_touch
    all_dat = get_s1s2_all_dat;
    ani_used = find(all_dat.has_deep_data_s1s2);

    % --- make the blank plots
    fh = figure('Position', [0 0 1000 500]);
    
    frac_by_type_and_area_ax = subplot('Position', [.1 .05 .2 .8]);
    touch_composition_ax = subplot('Position', [.4 .05 .2 .8]);

    % --- fraction of cells by area and type
    frac_all_cells_by_type_and_area = nan*zeros(length(all_dat.anims), 2, 3); % animal X s1/s2 X uSW/bSW/MW
    frac_touch_cells_by_type_and_area = nan*zeros(length(all_dat.anims), 2, 3); % animal X s1/s2 X uSW/bSW/MW
    d_s1 = 1; % most superficial of merged subvolumes
    d_s2 = 3; 
    for ai=1:length(ani_used)
        a = ani_used(ai);

        % figure out which neurons go where
%        s1i = find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s1_id_range{a}(2));
%        s2i = find(all_dat.anim_data{a}.ids > all_dat.s2_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s2_id_range{a}(2));
        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);
%        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_all', all_dat, a);

        % which date should we use? second non-compound date
        non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);

        d_s1_used = d_s1;
        d_s2_used = d_s2;
        % animal-specific stuff
        if (strcmp(all_dat.anims{a}, 'an014359'))
            d_s2_used = 1;
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
    end

    hold(frac_by_type_and_area_ax, 'on');
    hold(touch_composition_ax, 'on');
    cell_type = {'uSW','bSW','MW'};
    colors = {all_dat.usw_color, all_dat.bsw_color, all_dat.mw_color};
    for cti=1:3 % cell type loop
        v_s1 = squeeze(frac_all_cells_by_type_and_area(:,1,cti));
        v_s2 = squeeze(frac_all_cells_by_type_and_area(:,2,cti));
        [h pv] = ttest(v_s1,v_s2);
        pv = signrank(v_s1,v_s2);
        disp(sprintf('overall fraction for %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    cell_type{cti}, nanmean(v_s1), nanstd(v_s1), nanmean(v_s2), nanstd(v_s2), length(find(~isnan(v_s1))), pv));
        bar(frac_by_type_and_area_ax, cti+[0 4], [nanmean(v_s1) nanmean(v_s2)], .2, 'FaceColor', colors{cti}); 
        if (cti == 1)
            set(frac_by_type_and_area_ax, 'XTick', [2 6], 'XTickLabel', {'s1', 's2'});
            ylabel(frac_by_type_and_area_ax, 'Fraction of neurons');
        end
                    

        v_s1 = squeeze(frac_touch_cells_by_type_and_area(:,1,cti));
        v_s2 = squeeze(frac_touch_cells_by_type_and_area(:,2,cti));
        [h pv] = ttest(v_s1,v_s2);
        pv = signrank(v_s1,v_s2);
        disp(sprintf('touch only fraction for %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    cell_type{cti}, nanmean(v_s1), nanstd(v_s1), nanmean(v_s2), nanstd(v_s2), length(find(~isnan(v_s1))), pv));
        bar(touch_composition_ax, cti+[0 4], [nanmean(v_s1) nanmean(v_s2)], .2, 'FaceColor', colors{cti}); 
        if (cti == 1)
            set(touch_composition_ax, 'XTick', [2 6], 'XTickLabel', {'s1', 's2'});
            ylabel(touch_composition_ax, 'Fraction of touch neurons');
        end                    
    end

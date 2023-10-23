function ratio_mat= plot_proj_single_set(all_dat, vali, restrict_tag, ax, tstr, mode, ylab)
    
    s1_tag = 's1_all';
    %s1_tag = 's1_manual_restrict';
    %s2_tag = 's2_all';
    s2_tag = 's2_manual_restrict';




    if (nargin < 6) ; mode = 1; end % 1: fraction of projection 2: count relative expected by chance
    if (nargin < 7) ; ylab = '';  end
    cts = {'usw', 'bsw', 'mw'};

    ratio_mat = zeros(3,length(vali));
    for v=1:length(vali)
        ai=vali(v); 

        % restrict to desired area's id range
        area_nrni = get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, restrict_tag, all_dat, ai);

        % map dates onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{ai}.date_str));
        for d=1:length(all_dat.anim_data{ai}.date_str)
            if (length(all_dat.anim_data{ai}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}));
            end
        end
        di = min(find(all_dat.retro_dates{ai}));
        dat_dayi = find(d_num == di);

        % and now get the ratio of actual to expected fraction by type
        for c=1:length(cts) ; 
            ct = cts{c};
            cti = intersect(area_nrni,all_dat.types_by_idx{ai}.([ct '_by_day']){dat_dayi});
            N_ct(c) = length(cti);
            N_red_ct(c) = length(intersect(cti,all_dat.types_by_idx{ai}.is_red));
        end
        N_touch = sum(N_ct);
        N_all = length(area_nrni);
        N_red = length(all_dat.types_by_idx{ai}.is_red);

        resp_binary_mat = nan*zeros(4,length(all_dat.anim_data{ai}.ids));
        resp_binary_mat(1,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w1p_by_day{dat_dayi})) = 1;
        resp_binary_mat(2,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w1r_by_day{dat_dayi})) = 1;
        resp_binary_mat(3,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w2p_by_day{dat_dayi})) = 1;
        resp_binary_mat(4,intersect(area_nrni,all_dat.types_by_idx{ai}.touch_w1r_by_day{dat_dayi})) = 1;

        resp_dff_mat = nan*zeros(4,length(all_dat.anim_data{ai}.ids));
        resp_dff_mat(1,:) = all_dat.anim_mats{ai}.probRespW1p(dat_dayi,:);
        resp_dff_mat(2,:) = all_dat.anim_mats{ai}.probRespW1r(dat_dayi,:);
        resp_dff_mat(3,:) = all_dat.anim_mats{ai}.probRespW2p(dat_dayi,:);
        resp_dff_mat(4,:) = all_dat.anim_mats{ai}.probRespW2r(dat_dayi,:);

        mu_resp_per_cell = nanmean(resp_dff_mat.*resp_binary_mat);
%        mu_resp_per_cell = resp_dff_mat;    

        net_resp_proj = nansum(mu_resp_per_cell(all_dat.types_by_idx{ai}.is_red));
        for c=1:length(cts) ; 
            ct = cts{c};
            cti = intersect(area_nrni,all_dat.types_by_idx{ai}.([ct '_by_day']){dat_dayi});
            if (mode == 2) % number actual : number expected by chance (ratio)
                normd = N_red_ct(c)/((N_ct(c)/N_all)*(N_red/N_all)*N_all);
            elseif (mode == 1) % fraction of projection
                normd = N_red_ct(c)/N_red;
            elseif (mode == 3) % fraction of touchcells
                normd = N_red_ct(c)/N_touch; 
            elseif (mode == 4) % fraction of touchcells *of type*
                normd = N_red_ct(c)/N_ct(c);
            elseif (mode == 5) % raw cel lcount of type
                normd = N_red_ct(c);        
            elseif (mode == 6) % fraction of *activity* 
                pcti = intersect(cti, all_dat.types_by_idx{ai}.is_red);
                normd = nansum(mu_resp_per_cell(pcti))/net_resp_proj;
            elseif (mode == 7) % mean response, projecting
                pcti = intersect(cti, all_dat.types_by_idx{ai}.is_red);
                normd = nanmean(mu_resp_per_cell(pcti));
            elseif (mode == 8) % mean response, non-proj
                npcti = setdiff(cti, all_dat.types_by_idx{ai}.is_red);
                normd = nanmean(mu_resp_per_cell(npcti));
            end
        %    disp(sprintf('%s cell: %s expected: %0.2f real: %d ratio: %0.2f', all_dat.anims{ai}, ct, (N_ct(c)/N_all)*(N_red/N_all)*N_all, N_red_ct(c), normd));
            ratio_mat(c,v) = normd;
        end
        %disp( ' ' ); 
    end

    ms = 10;
    
    hold(ax, 'on');
    title(ax, tstr);

    xoffs = rand(size(ratio_mat,2),1)*0.2-0.1;
    xoffs = xoffs*0; % disable jitter
    
    if (1)
        x = 1;
        rectangle(ax, 'Position', [x-0.4 0 .8 nanmean(ratio_mat(1,:)) ], 'EdgeColor','None','FaceColor', all_dat.usw_color);
        plot(ax,0*ratio_mat(1,:)+x+xoffs, ratio_mat(1,:), 'o', 'Color', [1 1 1]*0.5, 'MarkerSize', ms, 'MarkerFaceColor', min(all_dat.usw_color+0.25, [1 1 1]), 'LineWidth', 2);

        x = 2;
        rectangle(ax, 'Position', [x-0.4 0 .8 nanmean(ratio_mat(2,:)) ], 'EdgeColor','None','FaceColor', all_dat.bsw_color);
        plot(ax,0*ratio_mat(1,:)+x+xoffs, ratio_mat(2,:), 'o', 'Color', [1 1 1]*0.5, 'MarkerSize', ms, 'MarkerFaceColor', min(all_dat.bsw_color+0.25, [1 1 1]), 'LineWidth', 2);

        x = 3;
        rectangle(ax, 'Position', [x-0.4 0 .8 nanmean(ratio_mat(3,:)) ], 'EdgeColor','None','FaceColor', all_dat.mw_color);
        plot(ax,0*ratio_mat(1,:)+x+xoffs, ratio_mat(3,:), 'o', 'Color', [1 1 1]*0.5, 'MarkerSize', ms, 'MarkerFaceColor', min(all_dat.mw_color+0.25, [1 1 1]), 'LineWidth', 2);

        [h pv_ab] = ttest2(ratio_mat(1,:), ratio_mat(2,:));
        [h pv_ac] = ttest2(ratio_mat(1,:), ratio_mat(3,:));
        [h pv_bc] = ttest2(ratio_mat(2,:), ratio_mat(3,:));

        title(ax, {tstr, sprintf('ab: %0.3f ac: %0.3f bc: %0.3f', pv_ab, pv_ac, pv_bc)});

        disp(sprintf('%s %s uSW: %0.2f/%0.2f (mu/sd) bSW: %0.2f/%0.2f MW: %0.2f/%0.2f pvals u v. b: %0.3f u v. m: %0.3f b v.m: %0.3f', ylab, tstr, ...
             nanmean(ratio_mat(1,:)), nanstd(ratio_mat(1,:)),nanmean(ratio_mat(2,:)), nanstd(ratio_mat(2,:)),nanmean(ratio_mat(3,:)), nanstd(ratio_mat(3,:)),  pv_ab, pv_ac, pv_bc));

    else
        plot(ax,0*ratio_mat(1,:)+1+xoffs, ratio_mat(1,:), 'o', 'Color', all_dat.usw_color, 'MarkerSize', ms);
        plot(ax,1, nanmean(ratio_mat(1,:)), 'o', 'Color', all_dat.usw_color, 'MarkerSize', ms+4, 'LineWidth',2);
        
        plot(ax,0*ratio_mat(2,:)+2+xoffs, ratio_mat(2,:), 'o', 'Color', all_dat.bsw_color, 'MarkerSize', ms);
        plot(ax,2, nanmean(ratio_mat(2,:)), 'o', 'Color', all_dat.bsw_color, 'MarkerSize', ms+4, 'LineWidth',2);

        plot(ax,0*ratio_mat(3,:)+3+xoffs, ratio_mat(3,:), 'o', 'Color', all_dat.mw_color, 'MarkerSize', ms);
        plot(ax,3, nanmean(ratio_mat(3,:)), 'o', 'Color', all_dat.mw_color, 'MarkerSize', ms+4, 'LineWidth',2);
    end    
    set(ax,'TickDir','out','FontSize',15);


    aa=axis(ax);
    axis(ax, [0 4 0 aa(4)]);
    ylabel(ax,ylab);


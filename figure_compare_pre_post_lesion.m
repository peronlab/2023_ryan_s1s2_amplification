%
% Figure for following panels:
%   - fraction of all cells by area that have specific type (uSW, bSW, MW)
%   - mean encoding score for touch cell (may want to take top 5% and mean of that, for example)
%   - ???
%
function figure_compare_pre_post_lesion
    all_dat = get_s1s2_all_dat;
    s1_lesion_ani = find(~strcmp('',all_dat.s1_lesion_date));
    s2_lesion_ani = find(~strcmp('',all_dat.s2_lesion_date));
%   sham_lesion_ani = find(all_dat.is_1daysham); % restrict to only 1day shams (This is more appropriate - what we do for everything else!)
    sham_lesion_ani = find(~strcmp('',all_dat.sham_lesion_date)); % use all shams

    % --- gather data pre/post
    [frac_all_cells_prepost_s1L presp_all_cells_prepost_s1L ddff_all_cells_prepost_s1L] = plot_touch_cells_pre_post_lesion (0, all_dat, s1_lesion_ani, all_dat.lesion_dates, 's2_manual_restrict', 'S1 lesion') ;   
    [frac_all_cells_prepost_s2L presp_all_cells_prepost_s2L ddff_all_cells_prepost_s2L] = plot_touch_cells_pre_post_lesion (600, all_dat, s2_lesion_ani, all_dat.lesion_dates, 's1_manual_restrict', 'S2 lesion');      
    [frac_all_cells_prepost_shL presp_all_cells_prepost_shL ddff_all_cells_prepost_shL] = plot_touch_cells_pre_post_lesion (1200, all_dat, sham_lesion_ani, all_dat.sham_dates, 's1_manual_restrict', 'Sham lesion');      

    delta_compare(frac_all_cells_prepost_s1L, frac_all_cells_prepost_shL, 's1-sh frac');
    delta_compare(presp_all_cells_prepost_s1L, presp_all_cells_prepost_shL, 's1-sh presp');
    delta_compare(ddff_all_cells_prepost_s1L, ddff_all_cells_prepost_shL, 's1-sh dfff');
    
    delta_compare(frac_all_cells_prepost_s2L, frac_all_cells_prepost_shL, 's2-sh frac');
    delta_compare(presp_all_cells_prepost_s2L, presp_all_cells_prepost_shL, 's2-sh presp');
    delta_compare(ddff_all_cells_prepost_s2L, ddff_all_cells_prepost_shL, 's2-sh dfff');
 
    plot_touch_cells_type_pre_post_lesion (0, all_dat, s1_lesion_ani, all_dat.lesion_dates, 's2_manual_restrict', 'S1 lesion');    
    plot_touch_cells_type_pre_post_lesion (600,all_dat, s2_lesion_ani, all_dat.lesion_dates, 's1_manual_restrict', 'S2 lesion');      
    plot_touch_cells_type_pre_post_lesion (1200, all_dat, sham_lesion_ani, all_dat.sham_dates, 's1_manual_restrict', 'Sham lesion');      

function delta_compare (prepost_1, prepost_2, tstr)
    delta_1 = (prepost_1(:,2)-prepost_1(:,1))./prepost_1(:,1);
    delta_2 = (prepost_2(:,2)-prepost_2(:,1))./prepost_2(:,1);

    disp(sprintf('comparing: %s pvalue: %0.3f', tstr, stattest_unpaired(delta_1,delta_2))); 

function plot_touch_cells_type_pre_post_lesion (x0, all_dat, ani_used, dates_used, restrict_tag, les_type)    
    % figure setup
    fh = figure ('Position', [x0 400 600 400]);
    raw_ax(1) = subplot('Position', [.1 .1 .2 .4]);
    delta_ax(1) = subplot('Position', [.1 .55 .2 .4]);
    raw_ax(2) = subplot('Position', [.4 .1 .2 .4]);
    delta_ax(2) = subplot('Position', [.4 .55 .2 .4]);
    raw_ax(3) = subplot('Position', [.7 .1 .2 .4]);
    delta_ax(3) = subplot('Position', [.7 .55 .2 .4]);

    disp(['----------------------------------' les_type '-----------------------------------------------']);

    frac_all_cells_prepost = nan*zeros(length(all_dat.anims), 2, 3); % animal X pre/post X usw/bsw/mw
    presp_all_cells_prepost = nan*zeros(length(all_dat.anims), 2, 3); % animal X pre/post X usw/bsw/mw
    ddff_all_cells_prepost = nan*zeros(length(all_dat.anims), 2, 3); % animal X pre/post X usw/bsw/mw
%    ani_used ; pause
    for ai=1:length(ani_used)
        a = ani_used(ai);
        ndays = length(all_dat.anim_data{a}.date_str);

        % restrict by ID 
        ni = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, restrict_tag, all_dat, a);
%        ni = find(all_dat.anim_data{a}.ids >= id_range{a}(1) & all_dat.anim_data{a}.ids <= id_range{a}(2));
        if (length(ni) == 0) ; disp([' **** PROBLEM **** ' all_dat.anims{a} ' has ZERO neurons ']) ; end
        
        % which date should we use?
        pre_di = find(dates_used{a} == -1);
        post_di = find(dates_used{a} == 1);

        % map onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{a}.date_str));
        for d=1:length(all_dat.anim_data{a}.date_str)
            if (length(all_dat.anim_data{a}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{a}, all_dat.anim_data{a}.date_str{d}));
            end
        end

        % --- fraction

        % for fraction, we compute *MEAN* if we have multiple pre/post dates
        frac_all_cells_prepost(a,1,:) = 0; 
        for d=1:length(pre_di)
            di = find(d_num == pre_di(d));
            usw_d = intersect(ni,all_dat.types_by_idx{a}.usw_by_day{di});
            bsw_d = intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{di});
            mw_d = intersect(ni,all_dat.types_by_idx{a}.mw_by_day{di});

            frac_all_cells_prepost(a,1,1) = frac_all_cells_prepost(a,1,1) + length(usw_d)/length(ni);
            frac_all_cells_prepost(a,1,2) = frac_all_cells_prepost(a,1,2) + length(bsw_d)/length(ni);
            frac_all_cells_prepost(a,1,3) = frac_all_cells_prepost(a,1,3) + length(mw_d)/length(ni);
        end
        frac_all_cells_prepost(a,1,1:3) = frac_all_cells_prepost(a,1,1:3)/length(pre_di);

        frac_all_cells_prepost(a,2,:) = 0; 
        for d=1:length(post_di)
            di = find(d_num == post_di(d));
            usw_d = intersect(ni,all_dat.types_by_idx{a}.usw_by_day{di});
            bsw_d = intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{di});
            mw_d = intersect(ni,all_dat.types_by_idx{a}.mw_by_day{di});

            frac_all_cells_prepost(a,2,1) = frac_all_cells_prepost(a,2,1) + length(usw_d)/length(ni);
            frac_all_cells_prepost(a,2,2) = frac_all_cells_prepost(a,2,2) + length(bsw_d)/length(ni);
            frac_all_cells_prepost(a,2,3) = frac_all_cells_prepost(a,2,3) + length(mw_d)/length(ni);
        end
        frac_all_cells_prepost(a,2,1:3) = frac_all_cells_prepost(a,2,1:3)/length(post_di);

        % --- response probability & mean dff -- we use neurons that were responsive pre AND/OR post
        uswi = [];
        bswi = [];
        mwi = [];
        w1pi = [];
        w1ri = [];
        w2pi = [];
        w2ri = [];

        for d=1:length(pre_di)
            di = find(d_num == pre_di(d));

            w1pi = union(w1pi, intersect(ni,all_dat.types_by_idx{a}.touch_w1p_by_day{di}));
            w1ri = union(w1ri, intersect(ni,all_dat.types_by_idx{a}.touch_w1r_by_day{di}));
            w2pi = union(w2pi, intersect(ni,all_dat.types_by_idx{a}.touch_w2p_by_day{di}));
            w2ri = union(w2ri, intersect(ni,all_dat.types_by_idx{a}.touch_w2r_by_day{di}));

            uswi = union(uswi, intersect(ni,all_dat.types_by_idx{a}.usw_by_day{di}));
            bswi = union(bswi, intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{di}));
            mwi = union(mwi, intersect(ni,all_dat.types_by_idx{a}.mw_by_day{di}));
        end

        for d=1:length(post_di)
            di = find(d_num == post_di(d));

            w1pi = union(w1pi, intersect(ni,all_dat.types_by_idx{a}.touch_w1p_by_day{di}));
            w1ri = union(w1ri, intersect(ni,all_dat.types_by_idx{a}.touch_w1r_by_day{di}));
            w2pi = union(w2pi, intersect(ni,all_dat.types_by_idx{a}.touch_w2p_by_day{di}));
            w2ri = union(w2ri, intersect(ni,all_dat.types_by_idx{a}.touch_w2r_by_day{di}));

            uswi = union(uswi, intersect(ni,all_dat.types_by_idx{a}.usw_by_day{di}));
            bswi = union(bswi, intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{di}));
            mwi = union(mwi, intersect(ni,all_dat.types_by_idx{a}.mw_by_day{di}));
        end
        whiski = {w1pi, w1ri, w2pi, w2ri};
        groupi = {uswi, bswi, mwi};
        
        % response probability -- use intersection of each type of touch with particular group
        % build a response matrix that nan's touch types that neuron does not respond to
        resp_mat_pre = nan*zeros(length(all_dat.anim_data{a}.ids), 4);
        dni = find(ismember(d_num, pre_di));
        resp_mat_pre = get_responsive_only_resp_mat (resp_mat_pre, all_dat, a, dni, whiski, 'prob');
        for g=1:length(groupi)
            presp_all_cells_prepost(a,1,g) = nanmean(nanmean(resp_mat_pre(intersect(groupi{g},unique([w1pi' w1ri' w2pi' w2ri'])),:)));
        end

        resp_mat_post = nan*zeros(length(all_dat.anim_data{a}.ids), 4);
        dni = find(ismember(d_num, post_di));
        resp_mat_post = get_responsive_only_resp_mat (resp_mat_post, all_dat, a, dni, whiski, 'prob');
        for g=1:length(groupi)
            presp_all_cells_prepost(a,2,g) = nanmean(nanmean(resp_mat_post(intersect(groupi{g},unique([w1pi' w1ri' w2pi' w2ri'])),:)));
        end

        resp_mat_pre = nan*resp_mat_pre;
        dni = find(ismember(d_num, pre_di));
        resp_mat_pre = get_responsive_only_resp_mat (resp_mat_pre, all_dat, a, dni, whiski, 'mean');
        for g=1:length(groupi)
            ddff_all_cells_prepost(a,1,g) = nanmean(nanmean(resp_mat_pre(intersect(groupi{g},unique([w1pi' w1ri' w2pi' w2ri'])),:)));
        end

        resp_mat_post = nan*resp_mat_post;
        dni = find(ismember(d_num, post_di));
        resp_mat_post = get_responsive_only_resp_mat (resp_mat_post, all_dat, a, dni, whiski, 'mean');
        for g=1:length(groupi)
            ddff_all_cells_prepost(a,2,g) = nanmean(nanmean(resp_mat_post(intersect(groupi{g},unique([w1pi' w1ri' w2pi' w2ri'])),:)));
        end
    end

    % --- bar plots comparing change 

    % normalized change in dff
    title(raw_ax(1), les_type);
    single_prepost_type_comparison(all_dat, [raw_ax(1) delta_ax(1)], frac_all_cells_prepost, 'touch fraction' , 'Type fraction', les_type);
    axis(raw_ax(1), [0 7 0 0.25]);
    axis(delta_ax(1), [0 4 -1 1.5]);
    if (strcmp(les_type, 'S1 lesion'))
        axis(raw_ax(1), [0 7 0 0.15]);
    end

    % normalized change in dff
    single_prepost_type_comparison(all_dat, [raw_ax(2) delta_ax(2)], presp_all_cells_prepost, 'response probability' , 'Mean response probability', les_type);
    axis(raw_ax(2), [0 7 0 0.4]);
    axis(delta_ax(2), [0 4 -1 1]);
    if (strcmp(les_type, 'S1 lesion'))
        axis(raw_ax(2), [0 7 0 0.3]);
    end

    % normalized change in dff
    single_prepost_type_comparison(all_dat, [raw_ax(3) delta_ax(3)], ddff_all_cells_prepost, 'ddff' , 'Mean \DeltaF/F', les_type);
    axis(raw_ax(3), [0 7 0 0.8]);
    axis(delta_ax(3), [0 4 -1 1]);
    if (strcmp(les_type, 'S1 lesion'))
        axis(raw_ax(3), [0 7 0 0.8]);
    end


function single_prepost_type_comparison(all_dat, ax, dat_mat, tstr, ylab, les_type)

    % gather data -- lesioned whisker if s1
    vv_pre = squeeze(dat_mat(:, 1, :));
    vv_post = squeeze(dat_mat(:, 2, :));
    ai = find(~isnan(vv_pre(:,1)));
    anis = {};
    for aa=1:length(ai)
        a = ai(aa);
        anis{aa} = all_dat.anims{a};
        pre_usw_vec(aa) = vv_pre(a,1);
        post_usw_vec(aa) = vv_post(a,1);
        pre_bsw_vec(aa) = vv_pre(a,2);
        post_bsw_vec(aa) = vv_post(a,2);
        pre_mw_vec(aa) = vv_pre(a,3);
        post_mw_vec(aa) = vv_post(a,3);
    end

    %axes(ax(1)); 
    hold(ax(1), 'on');
    hold(ax(2), 'on');
    ms=8;

    % aggregate
    cell_type = {'usw','bsw','mw'};
    plot(ax(2), [0 2*length(cell_type)+1], [0 0], 'k:');
    pre_vecs = {pre_usw_vec, pre_bsw_vec, pre_mw_vec};
    post_vecs = {post_usw_vec, post_bsw_vec, post_mw_vec};
    colrs = {all_dat.usw_color, all_dat.bsw_color, all_dat.mw_color};
    for cti=1:length(cell_type) 
        v_pre = pre_vecs{cti};
        v_post = post_vecs{cti};
        pv = stattest_paired(v_pre,v_post);
        disp(sprintf('%s for %s pre mean/sd: %0.3f/%0.3f post: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    tstr, cell_type{cti}, nanmean(v_pre), nanstd(v_pre), nanmean(v_post), nanstd(v_post), length(find(~isnan(v_pre))), pv)); 

        for a=1:length(v_post)
            if (~isnan(v_post(a)))
                plot(ax(1), (cti-1)*2 + [1 2], [v_pre(a) v_post(a)], '-', 'Color', min ([1 1 1], colrs{cti}+[0.25 0.25 0.25]));
            end
            plot(ax(1), (cti-1)*2+[1 2], [nanmean(v_pre) nanmean(v_post)], '-', 'Color', colrs{cti}, 'LineWidth',4);
        end                        
        ax(1).Title.String = [ax(1).Title.String ' ' sprintf('p: %0.3f', pv)];

        normd_delta = (v_post-v_pre)./v_pre;

%        bar(ax(2), cti, nanmean(normd_delta), 'FaceColor', colrs{cti});
        jit = 0.125*(rand(1,length(normd_delta))-0.5);
        jit=jit*0; % disable jitter
        plot(ax(2), 0*normd_delta+cti+jit, normd_delta, 'o', 'Color', [1 1 1]*0.25, 'MarkerFaceColor', [1 1 1]*0.75, 'MarkerSize', ms);
        plot(ax(2), [-0.25 0.25]+cti, [1 1]*nanmean(normd_delta), '-', 'Color', colrs{cti}, 'LineWidth', 5);
        ax(2).Title.String = [ax(2).Title.String ' ' sprintf('p: %0.3f', stattest_v_zero(normd_delta))];

        for a=1:length(v_post)
            if (~isnan(v_post(a)))
                %disp(sprintf('     %s pre: %0.3f post: %0.3f', anis{a}, v_pre(a), v_post(a)));
                disp(sprintf('     %s pre: %0.3f post: %0.3f delt: %0.3f', anis{a}, v_pre(a), v_post(a), normd_delta(a)));
            end
        end        
    end                        

    % cleanup
    ylabel(ax(1), ylab);
    set (ax(1),'TickDir','out','XTick',[]);
    set (ax(2),'TickDir','out','XTick',[]);

function [frac_all_cells_prepost presp_all_cells_prepost ddff_all_cells_prepost] = plot_touch_cells_pre_post_lesion (x0, all_dat, ani_used, dates_used, restrict_tag, les_type)    
    % figure setup
    fh = figure ('Position', [x0 0 600 400]);
    raw_ax(1) = subplot('Position', [.1 .1 .2 .4]);
    delta_ax(1) = subplot('Position', [.1 .55 .2 .4]);
    raw_ax(2) = subplot('Position', [.4 .1 .2 .4]);
    delta_ax(2) = subplot('Position', [.4 .55 .2 .4]);
    raw_ax(3) = subplot('Position', [.7 .1 .2 .4]);
    delta_ax(3) = subplot('Position', [.7 .55 .2 .4]);

    disp(['----------------------------------' les_type '-----------------------------------------------']);

    frac_all_cells_prepost = nan*zeros(length(all_dat.anims), 2, 3); % animal X pre/post X w1/w2/both
    presp_all_cells_prepost = nan*zeros(length(all_dat.anims), 2, 2); % animal X pre/post X w1/w2
    ddff_all_cells_prepost = nan*zeros(length(all_dat.anims), 2, 2); % animal X pre/post X w1/w2
    for ai=1:length(ani_used)
        a = ani_used(ai);
        ndays = length(all_dat.anim_data{a}.date_str);

        % restrict by ID 
        %ni = find(all_dat.anim_data{a}.ids >= id_range{a}(1) & all_dat.anim_data{a}.ids <= id_range{a}(2));
        ni = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, restrict_tag, all_dat, a);
        
        if (length(ni) == 0) ; disp([' **** PROBLEM **** ' all_dat.anims{a} ' has ZERO neurons ']) ; pause; end

        % which date should we use?
        pre_di = find(dates_used{a} == -1);
        post_di = find(dates_used{a} == 1);

        % map onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{a}.date_str));
        for d=1:length(all_dat.anim_data{a}.date_str)
            if (length(all_dat.anim_data{a}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{a}, all_dat.anim_data{a}.date_str{d}));
            end
        end

        % --- fraction

        % for fraction, we compute *MEAN* if we have multiple pre/post dates
        frac_all_cells_prepost(a,1,1:3) = 0; 
        for d=1:length(pre_di)
            di = find(d_num == pre_di(d));
            w1i_d = intersect(ni,union(all_dat.types_by_idx{a}.touch_w1p_by_day{di}, all_dat.types_by_idx{a}.touch_w1r_by_day{di}));
            w2i_d = intersect(ni,union(all_dat.types_by_idx{a}.touch_w2p_by_day{di}, all_dat.types_by_idx{a}.touch_w2r_by_day{di}));

            alli = union(w1i_d,w2i_d);
           
            w1ix_d = setdiff(w1i_d, w2i_d);
            w2ix_d = setdiff(w2i_d, w1i_d);
 
            frac_all_cells_prepost(a,1,1) = frac_all_cells_prepost(a,1,1) + length(w1ix_d)/length(ni);
            frac_all_cells_prepost(a,1,2) = frac_all_cells_prepost(a,1,2) + length(w2ix_d)/length(ni);
            frac_all_cells_prepost(a,1,3) = frac_all_cells_prepost(a,1,3) + length(alli)/length(ni);
        end
        frac_all_cells_prepost(a,1,1:3) = frac_all_cells_prepost(a,1,1:3)/length(pre_di);

        frac_all_cells_prepost(a,2,1:3) = 0; 
        for d=1:length(post_di)
            di = find(d_num == post_di(d));
            w1i_d = intersect(ni,union(all_dat.types_by_idx{a}.touch_w1p_by_day{di}, all_dat.types_by_idx{a}.touch_w1r_by_day{di}));
            w2i_d = intersect(ni,union(all_dat.types_by_idx{a}.touch_w2p_by_day{di}, all_dat.types_by_idx{a}.touch_w2r_by_day{di}));

            alli = union(w1i_d,w2i_d);

            w1ix_d = setdiff(w1i_d, w2i_d);
            w2ix_d = setdiff(w2i_d, w1i_d);

            frac_all_cells_prepost(a,2,1) = frac_all_cells_prepost(a,2,1) + length(w1ix_d)/length(ni);
            frac_all_cells_prepost(a,2,2) = frac_all_cells_prepost(a,2,2) + length(w2ix_d)/length(ni);
            frac_all_cells_prepost(a,2,3) = frac_all_cells_prepost(a,2,3) + length(alli)/length(ni);
        end
        frac_all_cells_prepost(a,2,1:3) = frac_all_cells_prepost(a,2,1:3)/length(post_di);

        % --- response probability -- we use neurons that were responsive pre AND/OR post
        presp_all_cells_prepost(a,1,1:2) = 0; 
        w1pi = [];
        w1ri = [];
        w2pi = [];
        w2ri = [];        
        for d=1:length(pre_di)
            di = find(d_num == pre_di(d));

            w1pi = union(w1pi, intersect(ni,all_dat.types_by_idx{a}.touch_w1p_by_day{di}));
            w1ri = union(w1ri, intersect(ni,all_dat.types_by_idx{a}.touch_w1r_by_day{di}));
            w2pi = union(w2pi, intersect(ni,all_dat.types_by_idx{a}.touch_w2p_by_day{di}));
            w2ri = union(w2ri, intersect(ni,all_dat.types_by_idx{a}.touch_w2r_by_day{di}));
        end
        for d=1:length(post_di)
            di = find(d_num == post_di(d));

            w1pi = union(w1pi, intersect(ni,all_dat.types_by_idx{a}.touch_w1p_by_day{di}));
            w1ri = union(w1ri, intersect(ni,all_dat.types_by_idx{a}.touch_w1r_by_day{di}));
            w2pi = union(w2pi, intersect(ni,all_dat.types_by_idx{a}.touch_w2p_by_day{di}));
            w2ri = union(w2ri, intersect(ni,all_dat.types_by_idx{a}.touch_w2r_by_day{di}));
        end
        whiski = {w1pi, w1ri, w2pi, w2ri};
        
        % build a response matrix that nan's touch types that neuron does not respond to
        resp_mat_pre = nan*zeros(length(all_dat.anim_data{a}.ids), 4);
        dni = find(ismember(d_num, pre_di));
        resp_mat_pre = get_responsive_only_resp_mat (resp_mat_pre, all_dat, a, dni, whiski, 'prob');
        presp_all_cells_prepost(a,1,1) = nanmean(nanmean(resp_mat_pre(union(w1pi,w1ri),[1 2])));
        presp_all_cells_prepost(a,1,2) = nanmean(nanmean(resp_mat_pre(union(w2pi,w2ri),[3 4])));

        resp_mat_post = nan*zeros(length(all_dat.anim_data{a}.ids), 4);
        dni = find(ismember(d_num, post_di));
        resp_mat_post = get_responsive_only_resp_mat (resp_mat_post, all_dat, a, dni, whiski, 'prob');
        presp_all_cells_prepost(a,2,1) = nanmean(nanmean(resp_mat_post(union(w1pi,w1ri),[1 2])));
        presp_all_cells_prepost(a,2,2) = nanmean(nanmean(resp_mat_post(union(w2pi,w2ri),[3 4])));

        % --- mean dff 
        resp_mat_pre = nan*resp_mat_pre;
        dni = find(ismember(d_num, pre_di));
        resp_mat_pre = get_responsive_only_resp_mat (resp_mat_pre, all_dat, a, dni, whiski, 'mean');
        ddff_all_cells_prepost(a,1,1) = nanmean(nanmean(resp_mat_pre(union(w1pi,w1ri),[1 2])));
        ddff_all_cells_prepost(a,1,2) = nanmean(nanmean(resp_mat_pre(union(w2pi,w2ri),[3 4])));

        resp_mat_post = nan*resp_mat_post;
        dni = find(ismember(d_num, post_di));
        resp_mat_post = get_responsive_only_resp_mat (resp_mat_post, all_dat, a, dni, whiski, 'mean');
        ddff_all_cells_prepost(a,2,1) = nanmean(nanmean(resp_mat_post(union(w1pi,w1ri),[1 2])));
        ddff_all_cells_prepost(a,2,2) = nanmean(nanmean(resp_mat_post(union(w2pi,w2ri),[3 4])));

    end

    % --- bar plots comparing change 

    % normalized change in dff
    title(raw_ax(1), les_type);
    single_prepost_comparison(all_dat, [raw_ax(1) delta_ax(1)], frac_all_cells_prepost, 'touch fraction' , 'Type fraction', les_type);
    axis(raw_ax(1), [0 3 0 0.3]);
    axis(delta_ax(1), [0 2 -1 1]);
    if (strcmp(les_type, 'S1 lesion'))
        axis(raw_ax(1), [0 7 0 0.2]);
        axis(delta_ax(1), [0 4 -1 1]);
    end

    % normalized change in dff
    single_prepost_comparison(all_dat, [raw_ax(2) delta_ax(2)], presp_all_cells_prepost, 'response probability' , 'Mean response probability', les_type);
    axis(raw_ax(2), [0 3 0 0.15]);
    axis(delta_ax(2), [0 2 -1 1]);
    if (strcmp(les_type, 'S1 lesion'))
        axis(raw_ax(2), [0 7 0 0.15]);
        axis(delta_ax(2), [0 4 -1 1]);
    end

    % normalized change in dff
    single_prepost_comparison(all_dat, [raw_ax(3) delta_ax(3)], ddff_all_cells_prepost, 'ddff' , 'Mean \DeltaF/F', les_type);
    axis(raw_ax(3), [0 3 0 0.25]);
    axis(delta_ax(3), [0 2 -1 1]);
    if (strcmp(les_type, 'S1 lesion'))
        axis(raw_ax(3), [0 7 0 0.25]);
        axis(delta_ax(3), [0 4 -1 1]);
    end

function resp_mat = get_responsive_only_resp_mat (resp_mat, all_dat, a, dni, whiski, mat_str)
        wstr = {'W1p','W1r','W2p','W2r'};

        for w=1:4
            field_str = [mat_str 'Resp' wstr{w}];
            if (length(dni) > 1)
                resp_mat(whiski{w},w) = nanmean(all_dat.anim_mats{a}.(field_str)(dni,whiski{w}));
            else
                resp_mat(whiski{w},w) = all_dat.anim_mats{a}.(field_str)(dni,whiski{w});
            end
        end


function single_prepost_comparison(all_dat, ax, dat_mat, tstr, ylab, les_type)

    % gather data -- lesioned whisker if s1
    whi = all_dat.s1_lesion_whisker;
    nwhi = whi;
    nwhi(find(whi == 2)) = 1;
    nwhi(find(whi == 1)) = 2;
    vv_pre = squeeze(dat_mat(:, 1, 1:2));
    vv_post = squeeze(dat_mat(:, 2, 1:2));
    ai = find(~isnan(vv_pre(:,1)));
    anis = {};
    for aa=1:length(ai)
        a = ai(aa);
        anis{aa} = all_dat.anims{a};
        if (strcmp(tstr, 'touch fraction')) % here we sum since we want 
            pre_vec(aa) = dat_mat(a,1,3);
            post_vec(aa) = dat_mat(a,2,3);
        else
            pre_vec(aa) = nanmean(vv_pre(a,:));
            post_vec(aa) = nanmean(vv_post(a,:));
        end

        % separate lesioned/nonlesioned
        if (ismember(whi(a),[1 2]))
            pre_lesioned_vec(aa) = vv_pre(a,whi(a));
            pre_unlesioned_vec(aa) = vv_pre(a,nwhi(a));
            post_lesioned_vec(aa) = vv_post(a,whi(a));
            post_unlesioned_vec(aa) = vv_post(a,nwhi(a));            
        else
            pre_lesioned_vec(aa) = nan;
            pre_unlesioned_vec(aa) = nan;
            post_lesioned_vec(aa) = nan;
            post_unlesioned_vec(aa) = nan;
        end
    end

    %axes(ax(1)); 
    hold(ax(1), 'on');
    hold(ax(2), 'on');

    % aggregate
    mu_pre = nanmean(pre_vec);
    sd_pre = nanstd(pre_vec);
    mu_post = nanmean(post_vec);
    sd_post = nanstd(post_vec);
    disp(sprintf('%s for both whiskers pre: %0.3f/%0.3f mu/SE post: %0.3f/%0.3f; vs: %0.3f n=%d', tstr, mu_pre, sd_pre, mu_post, sd_post, stattest_paired(pre_vec,post_vec),length(find(~isnan(pre_vec)))));
    ax(1).Title.String = [ax(1).Title.String ' ' sprintf('p: %0.3f', stattest_paired(pre_vec,post_vec))];
    delta = (post_vec-pre_vec);
    disp(sprintf('              delta: %0.3f/%0.3f; vs 0: %0.3f', nanmean(delta), nanstd(delta), stattest_v_zero(delta)));
    normd_delta = (post_vec-pre_vec)./pre_vec;
    disp(sprintf('              normd delta: %0.3f/%0.3f; vs 0: %0.3f', nanmean(normd_delta), nanstd(normd_delta), stattest_v_zero(normd_delta)));
    
    for a=1:length(post_vec)
        if (~isnan(post_vec(a)))
            disp(sprintf('     %s pre: %0.3f post: %0.3f delt: %0.3f', anis{a}, pre_vec(a), post_vec(a), normd_delta(a)));
      %      disp(sprintf('     %s pre: %0.3f post: %0.3f', anis{a}, pre_vec(a), post_vec(a)));
            plot(ax(1), [1 2], [pre_vec(a) post_vec(a)], '-', 'Color', [ 1 1 1]*.5);
            plot(ax(1), [1 2], [nanmean(pre_vec) nanmean(post_vec)], 'k-', 'LineWidth', 4);
%            plot(ax(1), (cti-1)*2 + [1 2], [v_pre(a) v_post(a)], '-', 'Color', min ([1 1 1], colrs{cti}+[0.25 0.25 0.25]));
        end
    end
    ms=8;
%    bar(ax(2), 1, nanmean(normd_delta), 'FaceColor', [0 0 0]);
    plot(ax(2), 1+[-0.25 0.25], [1 1]*nanmean(normd_delta), '-', 'Color', [0 0 0], 'LineWidth', 5);
    jit = 0.125*(rand(1,length(normd_delta))-0.5);
    jit=jit*0; % disable jitter
    
    plot(ax(2), 0*normd_delta+1+jit, normd_delta, 'o', 'Color', [1 1 1]*0.25, 'MarkerFaceColor', [1 1 1]*0.75, 'MarkerSize', ms);
    title(ax(2), sprintf('p: %0.3f', stattest_v_zero(normd_delta)));

    if (strcmp(les_type, 'S1 lesion')) % we will plot lesioned and unlesioned in this case
        cell_type = {'lesioned','unlesioned'};
        pre_vecs = {pre_lesioned_vec, pre_unlesioned_vec};
        post_vecs = {post_lesioned_vec, post_unlesioned_vec};
        normd_delta = {};
        for cti=1:length(cell_type) 
            v_pre = pre_vecs{cti};
            v_post = post_vecs{cti};
            pv = stattest_paired(v_pre,v_post);
            disp(sprintf('%s for %s pre mean/sd: %0.3f/%0.3f post: %0.3f/%0.3f n=%d pval: %0.3f', ...
                        tstr, cell_type{cti}, nanmean(v_pre), nanstd(v_pre), nanmean(v_post), nanstd(v_post), length(find(~isnan(v_pre))), pv)); 
            ax(1).Title.String = [ax(1).Title.String ' ' sprintf('p: %0.3f', pv)];

            if (cti == 1) ; colr = [1 0 0] ; else ; colr = [1 1 1]*0.5 ; end

            for a=1:length(v_post)
                if (~isnan(v_post(a)))
                    plot(ax(1), cti*2 + [1 2], [v_pre(a) v_post(a)], '-', 'Color', min ([1 1 1], colr+[0.25 0.25 0.25]));
                    plot(ax(1), cti*2 + [1 2], [nanmean(v_pre) nanmean(v_post)], '-', 'Color', colr, 'LineWidth', 4);
                    
                end
            end                        
            
            A= (post_lesioned_vec-pre_lesioned_vec)./(pre_lesioned_vec);
            B= (post_unlesioned_vec-pre_unlesioned_vec)./(pre_unlesioned_vec);
            pv= stattest_paired(A,B);
            disp(sprintf('%s change lesioned versus unlesioned n=%d pval: %0.3f', ...
                        tstr, length(find(~isnan(v_pre))), pv)); 
            
            normd_delta{cti} = (v_post-v_pre)./v_pre;
            jit = 0.125*(rand(1,length(normd_delta{cti}))-0.5);
            jit=jit*0; % disable jitter
%            bar(ax(2), cti+1, nanmean(normd_delta{cti}), 'FaceColor', colr);
            plot(ax(2), cti+1+[-0.25 0.25], [1 1]*nanmean(normd_delta{cti}), '-', 'Color', colr, 'LineWidth', 5);
            plot(ax(2), 0*normd_delta{cti}+cti+1+jit, normd_delta{cti}, 'o', 'Color', [1 1 1]*0.25, 'MarkerFaceColor', [1 1 1]*0.75, 'MarkerSize', ms);
            ax(2).Title.String = [ax(2).Title.String ' ' sprintf('p: %0.3f', stattest_v_zero(normd_delta{cti}))];
        end                        

        ax(2).Title.String = [ax(2).Title.String ' ' sprintf('p L v nL: %0.3f', stattest_paired(normd_delta{1},normd_delta{2}))];
        plot(ax(2), [0 7], [0 0], 'k:');
    else
        plot(ax(2), [0 3], [0 0], 'k:');
    end

    % cleanup
    ylabel(ax(1), ylab);
    set (ax(1),'TickDir','out','XTick',[]);
    set (ax(2),'TickDir','out','XTick',[]);


function pv = stattest_unpaired(vec1, vec2)
    if (length(find(isnan(vec1))) == length(vec1))
        pv = -1; % meaningless
    else
        pv = ranksum(vec1, vec2);
    end
    [h pv] = ttest2(vec1, vec2);

function pv = stattest_paired(vec1, vec2)
    if (length(find(isnan(vec1))) == length(vec1))
        pv = -1; % meaningless
    else
        pv = signrank(vec1, vec2);
    end
    [h pv] = ttest(vec1, vec2);

function pv = stattest_v_zero(vec)
    if (length(find(isnan(vec))) == length(vec))
        pv = -1; % meaningless
    else
        pv = signrank(vec);
    end
    [h pv] = ttest(vec);

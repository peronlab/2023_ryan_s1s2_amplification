% Script for stats comparing proj and not. For figure 3 non correlation
% stats (mean)
% calls function: plot_proj_single_set

%% setup, then call fnctn -- changed the function to output ratio mat for each type, so we can compare across proj and non proj. 
    fh = figure('Position',[0 0 1500 900]);
    s1_to_s2_normd_ax = subplot('Position', [.1 .1 .15 .25]);
    s2_to_s1_normd_ax = subplot('Position', [.3 .1 .15 .25]);
    s1_to_s2_abs_ax = subplot('Position', [.1 .4 .15 .25]);
    s2_to_s1_abs_ax = subplot('Position', [.3 .4 .15 .25]);
    s1_to_s2_frac_ax = subplot('Position', [.1 .7 .15 .25]);
    s2_to_s1_frac_ax = subplot('Position', [.3 .7 .15 .25]);
    s1_to_s2_frac_of_subtype_ax = subplot('Position', [.55 .7 .15 .25]);
    s2_to_s1_frac_of_subtype_ax = subplot('Position', [.75 .7 .15 .25]);

    s1_to_s2_mean_resp_proj_ax = subplot('Position', [.55 .4 .15 .25]);
    s2_to_s1_mean_resp_proj_ax = subplot('Position', [.75 .4 .15 .25]);
    s1_to_s2_mean_resp_nonproj_ax = subplot('Position', [.55 .1 .15 .25]);
    s2_to_s1_mean_resp_nonproj_ax = subplot('Position', [.75 .1 .15 .25]);



    all_dat = get_s1s2_all_dat;

    s1_tag = 's1_all';
    %s1_tag = 's1_manual_restrict';
    %s2_tag = 's2_all';
    s2_tag = 's2_manual_restrict';

    vali = find(all_dat.retrograde_label_type ==  1); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    s1p_nonproj= plot_proj_single_set(all_dat, vali, s2_tag, s1_to_s2_mean_resp_proj_ax, 's1p neurons in s2', 7, 'Nean response to touch, proj (dFF)');

    vali = find(all_dat.retrograde_label_type ==  2); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    s2p_nonproj=plot_proj_single_set(all_dat, vali, s1_tag, s2_to_s1_mean_resp_proj_ax, 's2p neurons in s1', 7);

    aa_12 = axis(s1_to_s2_mean_resp_proj_ax);
    aa_21 = axis(s2_to_s1_mean_resp_proj_ax);
    axis(s1_to_s2_mean_resp_proj_ax, [aa_12(1:3) max(aa_12(4),aa_21(4))]);
    axis(s2_to_s1_mean_resp_proj_ax, [aa_21(1:3) max(aa_12(4),aa_21(4))]);
    M = max(aa_12(4),aa_21(4));
    
    vali = find(all_dat.retrograde_label_type ==  1); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    s1p_proj=plot_proj_single_set(all_dat, vali, s2_tag, s1_to_s2_mean_resp_nonproj_ax, 's1p neurons in s2', 8, 'Nean response to touch, nonproj (dFF)');

    vali = find(all_dat.retrograde_label_type ==  2); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    s2p_proj=plot_proj_single_set(all_dat, vali, s1_tag, s2_to_s1_mean_resp_nonproj_ax, 's2p neurons in s1', 8);

    aa_12 = axis(s1_to_s2_mean_resp_nonproj_ax);
    aa_21 = axis(s2_to_s1_mean_resp_nonproj_ax);
    axis(s1_to_s2_mean_resp_nonproj_ax, [aa_12(1:3) M]);
    axis(s2_to_s1_mean_resp_nonproj_ax, [aa_21(1:3) M]);
    %% STATS! 
    
    for aa=1:3
        A= s1p_nonproj(aa,:);
        B= s1p_proj(aa,:);
        [h,pv]= ttest(A,B);
        p_val(aa)= pv;
    end 
    
 
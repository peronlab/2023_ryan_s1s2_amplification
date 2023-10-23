function fig_kinematics_change
    all_dat = get_s1s2_all_dat;
    root_dir= '/Volumes/volume_imaging/';
    %% Example mouse plots

    % example mouse/day
    %ex_ai = find(strcmp(all_dat.anims, 'an018920')); % must be w/ whiskvid
    %ex_di = min(find(~isnan(all_dat.anim_dat(ex_ai).frac_corr_t))); % use pre-lesion day

    % setup the plots

    % pre/post, over several days, trial distro of mean/max amplitude, mean setpoint pretouch ; also look @ effect of sham

    % touch count pre/post

    % net dK pre/post *at specific positions* (distal and prox) on touch trials
    
    %% aggregate kinematic changes across mice
   s1_lesion_ani = find(~strcmp('',all_dat.s1_lesion_date));
   %s1_lesion_ani=s1_lesion_ani(1:2);
   s2_lesion_ani = find(~strcmp('',all_dat.s2_lesion_date));
   s1_list=  all_dat.anims(s1_lesion_ani);
   s2_list= all_dat.anims(s2_lesion_ani); 


    plot_aggregate_kinematic_changes (root_dir, all_dat, s1_lesion_ani, s1_list, all_dat.lesion_dates, 'S1 lesion' )

    plot_aggregate_kinematic_changes (root_dir, all_dat, s2_lesion_ani,s2_list, all_dat.lesion_dates, 'S2 lesion' )

    %% kinematic change as function of lesion size (LR)
    
function plot_aggregate_kinematic_changes (root_dir, all_dat, all_ai,  ani_list,dates_used, type_str)
    % setup the plots
    fh = figure('Position', [0 0 1500 600]);
    dx = 1/6;
    sp(1) = subplot('Position', [.05 .55 .12 .3]); 
    sp(2) = subplot('Position', [.05+dx .55 .12 .3]);
    sp(3) = subplot('Position', [.05+2*dx .55 .12 .3]);
    sp(4) = subplot('Position', [.05+3*dx .55 .12 .3]);
    sp(5) = subplot('Position', [.05+4*dx .55 .12 .3]);
    sp(6) = subplot('Position', [.05+5*dx .55 .12 .3]);
    sp(7) = subplot('Position', [.05 .075 .12 .3]);
    sp(8) = subplot('Position', [.05+dx .075 .12 .3]);
    sp(9) = subplot('Position', [.05+2*dx .075 .12 .3]);
    sp(10) = subplot('Position', [.05+3*dx .075 .12 .3]);
    sp(11) = subplot('Position', [.05+4*dx .075 .12 .3]);
    sp(12) = subplot('Position', [.05+5*dx .075 .12 .3]);

    for s=1:length(sp)
        set(sp(s), 'TickDir', 'out', 'FontSize', 15);
        hold(sp(s),'on');
    end
 
    
  for v=1:length(all_ai)
       
    net_dk_prelick = nan*zeros(length(all_ai), 2);
    q99_dk_prelick = nan*zeros(length(all_ai), 2);
    q99_amp_prelick = nan*zeros(length(all_ai), 2);
    q99_sp_prelick = nan*zeros(length(all_ai), 2);
    q99_vel_prelick = nan*zeros(length(all_ai), 2);
    mean_theta= nan*zeros(length(all_ai),2);
      
    ai=all_ai(v);
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
        
        anim= ani_list{v};
        cd(sprintf('%s/%s/session_neuropilone/',root_dir, anim)); %go to animal directory 
        date_list_pre=  all_dat.anim_data{ai}.date_str{day_pre};
        date_list_post= all_dat.anim_data{ai}.date_str{day_post};
        load (sprintf('%s_2022_%s_sess.mat', anim, date_list_pre));
        wdat_pre= s; 
        load (sprintf('%s_2022_%s_sess.mat', anim, date_list_post));
        wdat_post=s; 
        
        stats_pre= get_session_beh_stats(wdat_pre); 
        stats_post= get_session_beh_stats(wdat_post);
        
        num_touch_pre (v,1)= length(wdat_pre.whiskerBarContactESA.esa{1}.eventTimes)/length(wdat_pre.trial);%avg number of touches per trial
        num_touch_post (v,1)= length(wdat_post.whiskerBarContactESA.esa{1}.eventTimes)/length(wdat_post.trial);%avg number of touches per trial
        num_touch_pre (v,2)= length(wdat_pre.whiskerBarContactESA.esa{2}.eventTimes)/length(wdat_pre.trial);%avg number of touches per trial
        num_touch_post (v,2)= length(wdat_post.whiskerBarContactESA.esa{2}.eventTimes)/length(wdat_post.trial);%avg number of touches per trial
        
         
      
        %net_dk_prelick(a,1) = nanmedian([ ld.anim_dat(ai).trial_dat{di(d)}.net_kappa_pre_first_lick]);
        q99_dk_pre(v,1) = nanmedian( [stats_pre.trial_dat.q99_kappa_pre_first_lick]);
        q99_dk_post(v,1) = nanmedian( [stats_post.trial_dat.q99_kappa_pre_first_lick]);
        
        mean_theta_pre(v,1)= nanmedian( [stats_pre.trial_dat.mean_theta_at_touch]);
        mean_theta_post (v,1)= nanmedian( [stats_post.trial_dat.mean_theta_at_touch]);
        q99_amp_pre(v,1) = nanmedian( [stats_pre.trial_dat.q99_amplitude_pre_first_lick]);
        q99_amp_post(v,1) = nanmedian( [stats_post.trial_dat.q99_amplitude_pre_first_lick]);

        %q99_sp_prelick(a,1) = nanmedian( [ld.anim_dat(ai).trial_dat{di(d)}.q99_setpoint_pre_first_lick]);
        q99_vel_pre(v,1) = nanmedian([stats_pre.trial_dat.q99_velocity_pre_first_lick]);
        q99_vel_post(v,1) = nanmedian([stats_post.trial_dat.q99_velocity_pre_first_lick]);

        
        
        
  end   
    %%
    % plot - touch count
    
    for a=1:length(all_ai)
        ai = all_ai(a);
        plot(sp(1), [1 2], [num_touch_pre(a,1), num_touch_post(a,1)], '-','LineWidth',2, 'Color', [0.5 0.5 0.5]); hold on;
        %plot(sp(1), 2, num_touch_post(:,1), 'o-','LineWidth',2, 'Color', [0.5 0.5 0.5]);
        plot(sp(1),[1 2], [mean(num_touch_pre(:,1)), mean(num_touch_post(:,1))], '-o', 'LineWidth', 5, 'Color', [0.5 0.8 0.5]);
        plot(sp(2), [1 2], [num_touch_pre(a,2),num_touch_post(a,2)],  '-','LineWidth',2, 'Color', [0.25 0.35 0.25]); 
        %plot(sp(2), 2, num_touch_post(:,2), 'o-','LineWidth',2, 'Color', [0.25 0.35 0.25]);
        plot(sp(2),[1 2], [mean(num_touch_pre(:,2)), mean(num_touch_post(:,2))], '-o', 'LineWidth', 5, 'Color', [0.25 0.8 0.25]);
        
    end

   
    axis(sp(1), [0 3 0 14]);
    axis(sp(2), [0 3 0 14]);
    
    pv = -1;
    a = num_touch_pre(:,1);
    b = num_touch_post(:,1);
    if (length(find(~isnan(a+b))) > 1)
        [h pv] = ttest(a,b);
    end

    %disp(sprintf('num touch per trial %s %s pre: %0.2f/%0.2f post: %0.2f %0.2f', type_str, ev_descr_str, nanmean(a), nanstd(a), nanmean(b), nanstd(b)));

    %xlabel(sp(1), ['Days since ' ev_descr_str]);
    ylabel(sp(1), ['Number of touches c2 ']);
    title(sp(1),sprintf('%s; p_0_v_-_1: %0.3f ', type_str, pv));
   
    pv = -1;
    a = num_touch_pre(:,2);
    b = num_touch_post(:,2);
    if (length(find(~isnan(a+b))) > 1)
        [h pv] = ttest(a,b);
    end

    ylabel(sp(2), ['Number of touches c3']);
    title(sp(2),sprintf('%s; p_0_v_-_1: %0.3f ', type_str, pv));

%%
    % plot - max dk
    %xoffs = linspace(-0.1, 0.1, length(all_ai));
    for a=1:length(all_ai)
        ai = all_ai(a);
        plot(sp(3),[1 2], [mean(q99_dk_pre(:,1)), mean(q99_dk_post(:,1))], '-o', 'LineWidth', 5, 'Color', [0.5 0.8 0.5]);
        plot(sp(3), [1 2], [q99_dk_pre(a,1),q99_dk_post(a,1)],'-', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
        %plot(sp(3), 2, q99_dk_post(:,1), '--', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
    end

    axis(sp(3), [0 3 0 0.002]);
    
    pv = -1;
    a = q99_dk_pre(:,1);
    b = q99_dk_post(:,1);
    if (length(find(~isnan(a+b))) > 1)
        [h pv] = ttest(a,b);
    end

    %xlabel(sp(3), ['Days since ' ev_descr_str]);
    %isp(sprintf('99%% dK %s %s pre: %0.4f/%0.4f post: %0.4f %0.4f', type_str, ev_descr_str, nanmean(a), nanstd(a), nanmean(b), nanstd(b)));
    ylabel(sp(3), ['99% dK']);
    sp(3).YAxis.Exponent = 0;
    ytickformat(sp(3), '%03.3f');
    %ylim([0 0.001])
    title(sp(3),sprintf('%s; p_0_v_-_1: %0.3f ', type_str, pv));
%%
    % plot - mean theta
    %xoffs = linspace(-0.1, 0.1, length(all_ai));
    for a=1:length(all_ai)
        ai = all_ai(a);
        plot(sp(4),[1 2], [mean(mean_theta_pre(:,1)), mean(mean_theta_post(:,1))], '-o', 'LineWidth', 5, 'Color', [0.5 0.8 0.5]);
        plot(sp(4), [1 2], [mean_theta_pre(a,1),mean_theta_post(a,1)], '-', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
        %plot(sp(4), 2, mean_theta_post(:,1), '--', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
    end

   
    %axis(sp(3), [0 3 0 0.005]);
    axis(sp(4), [0 3 -20 0]);

    pv = -1;
    a = mean_theta_pre(:,1);
    b = mean_theta_post(:,1);
    if (length(find(~isnan(a+b))) > 1)
        [h pv] = ttest(a,b);
    end

    %xlabel(sp(3), ['Days since ' ev_descr_str]);
    %isp(sprintf('99%% dK %s %s pre: %0.4f/%0.4f post: %0.4f %0.4f', type_str, ev_descr_str, nanmean(a), nanstd(a), nanmean(b), nanstd(b)));
    ylabel(sp(4), ['mean theta']);
    sp(4).YAxis.Exponent = 0;
    ytickformat(sp(4), '%03.3f');
    title(sp(4),sprintf('%s; p_0_v_-_1: %0.3f ', type_str, pv));

% %%
%%
    % plot - velocity max
    %xoffs = linspace(-0.1, 0.1, length(all_ai));
    for a=1:length(all_ai)
        ai = all_ai(a);
        plot(sp(5),[1 2], [mean(q99_vel_pre(:,1)), mean(q99_vel_post(:,1))], '-o', 'LineWidth', 5, 'Color', [0.5 0.8 0.5]);
        plot(sp(5), [1 2], [q99_vel_pre(a,1), q99_vel_post(a,1)], '-', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
        %plot(sp(5), 2, q99_vel_post(:,1), '--', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
    end

   
    %axis(sp(3), [0 3 0 0.005]);
    axis(sp(5), [0 3 0 1000]);

    pv = -1;
    a = q99_vel_pre(:,1);
    b = q99_vel_post(:,1);
    if (length(find(~isnan(a+b))) > 1)
        [h pv] = ttest(a,b);
    end

    %xlabel(sp(3), ['Days since ' ev_descr_str]);
    %isp(sprintf('99%% dK %s %s pre: %0.4f/%0.4f post: %0.4f %0.4f', type_str, ev_descr_str, nanmean(a), nanstd(a), nanmean(b), nanstd(b)));
    ylabel(sp(5), ['99% velocity']);
    sp(5).YAxis.Exponent = 0;
    ytickformat(sp(5), '%03.3f');
    title(sp(5),sprintf('%s; p_0_v_-_1: %0.3f ', type_str, pv));
    %% 
        % plot - amplitude max
    %xoffs = linspace(-0.1, 0.1, length(all_ai));
    for a=1:length(all_ai)
        ai = all_ai(a);
        plot(sp(6),[1 2], [mean(q99_amp_pre(:,1)), mean(q99_amp_post(:,1))], '-o', 'LineWidth', 5, 'Color', [0.5 0.8 0.5]);
        plot(sp(6), [1 2], [q99_amp_pre(a,1), q99_amp_post(a,1)], '-', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
        %plot(sp(6), 2, q99_amp_post(:,1), '--', 'LineWidth',2, 'Color', [0.5 0.5 0.5]);
    end

   
    %axis(sp(3), [0 3 0 0.005]);
    axis(sp(6), [0 3 0 10]);

    pv = -1;
    a = q99_amp_pre(:,1);
    b = q99_amp_post(:,1);
    if (length(find(~isnan(a+b))) > 1)
        [h pv] = ttest(a,b);
    end

    %xlabel(sp(3), ['Days since ' ev_descr_str]);
    %isp(sprintf('99%% dK %s %s pre: %0.4f/%0.4f post: %0.4f %0.4f', type_str, ev_descr_str, nanmean(a), nanstd(a), nanmean(b), nanstd(b)));
    ylabel(sp(6), ['99% amplitdue']);
    sp(6).YAxis.Exponent = 0;
    ytickformat(sp(6), '%03.3f');
    title(sp(6),sprintf('%s; p_0_v_-_1: %0.3f ', type_str, pv));
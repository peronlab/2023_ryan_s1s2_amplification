%
% wrapper for plot_3d_s1s2
%
function fh= figure_s1s2_volume_maps
    fh=figure('Position',[0 0 300 600]);

    s1_w1p_ax = axes('Position', [0.1 0.55 0.4 0.4]);
    s2_w1p_ax = axes('Position', [0.55 0.55 0.4 0.4]);
    s1_w2p_ax = axes('Position', [0.1 0.1 0.4 0.4]);
    s2_w2p_ax = axes('Position', [0.55 0.1 0.4 0.4]);

    all_dat = get_s1s2_all_dat;
    %all_dat=all_dat;
    % parameters for plot_3d
    ani_str = 'an017510'%'an018920' % F4'an018927'; %an016650, F2 %S1 example an: an017510 day 5. S2 example an: 16650 6
    day_idx= 5;%5; %F2
    %day_idx_pre = %1 ;%2; %lesion, f4
    %day_idx_post= %3 ;%5; %lesion, f4
    %day_idx_sham= 3; %lesion, f4
    max_presp = 0.5;
    ball_size_range = [1 10];
    
%     %% section for pre/post lesion plotting figure 4/5
%     
%      plot_3d_s1s2(all_dat, ani_str, day_idx_pre, s1_w1p_ax, 'usw', 'w1p', 's2_manual_restrict', 1, max_presp, ball_size_range);
%      plot_3d_s1s2(all_dat, ani_str, day_idx_pre, s1_w1p_ax, 'bsw', 'w1p', 's2_manual_restrict', 0, max_presp, ball_size_range);
%      plot_3d_s1s2(all_dat, ani_str, day_idx_pre, s1_w1p_ax, 'mw', 'w1p', 's2_manual_restrict', 0, max_presp, ball_size_range);
%     
%      plot_3d_s1s2(all_dat, ani_str, day_idx_pre, s2_w1p_ax, 'usw', 'w2p', 's2_manual_restrict', 1, max_presp, ball_size_range);
%      plot_3d_s1s2(all_dat, ani_str, day_idx_pre, s2_w1p_ax, 'bsw', 'w2p', 's2_manual_restrict', 0, max_presp, ball_size_range);
%      plot_3d_s1s2(all_dat, ani_str, day_idx_pre, s2_w1p_ax, 'mw', 'w2p', 's2_manual_restrict', 0, max_presp, ball_size_range);
% %    
%        plot_3d_s1s2(all_dat, ani_str, day_idx_post, s1_w2p_ax, 'usw', 'w1p', 's2_manual_restrict', 1, max_presp, ball_size_range);
%        plot_3d_s1s2(all_dat, ani_str, day_idx_post, s1_w2p_ax, 'bsw', 'w1p', 's2_manual_restrict', 0, max_presp, ball_size_range);
%        plot_3d_s1s2(all_dat, ani_str, day_idx_post, s1_w2p_ax, 'mw', 'w1p', 's2_manual_restrict', 0, max_presp, ball_size_range);
% % %     
%        plot_3d_s1s2(all_dat, ani_str, day_idx_post, s2_w2p_ax, 'usw', 'w2p', 's2_manual_restrict',1, max_presp, ball_size_range);
%        plot_3d_s1s2(all_dat, ani_str, day_idx_post, s2_w2p_ax, 'bsw', 'w2p', 's2_manual_restrict', 0, max_presp, ball_size_range);
%        plot_3d_s1s2(all_dat, ani_str, day_idx_post, s2_w2p_ax, 'mw', 'w2p', 's2_manual_restrict', 0, max_presp, ball_size_range);
% % % %    
     %% Section for normal plotting of cells for figure 2
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'usw', 'w1p', 's1_manual_restrict', 1, max_presp, ball_size_range);
    plot_3d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'bsw', 'w1p', 's1_manual_restrict', 1, max_presp, ball_size_range);
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s1_w1p_ax, 'mw', 'w1p', 's1_manual_restrict', 0, max_presp, ball_size_range);

    %plot_3d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'usw', 'w1p', 's2_manual_restrict', 0, max_presp, ball_size_range);
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'bsw', 'w1p', 's2_manual_restrict', 0, max_presp, ball_size_range);
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s2_w1p_ax, 'mw', 'w1p', 's2_manual_restrict', 1, max_presp, ball_size_range);
    
   %plot_3d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'usw', 'w2p', 's1_manual_restrict', 1, max_presp, ball_size_range);
    plot_3d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'bsw', 'w2p', 's1_manual_restrict',0, max_presp, ball_size_range);
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s1_w2p_ax, 'mw', 'w2p', 's1_manual_restrict', 0, max_presp, ball_size_range);

   %^ plot_3d_s1s2(all_dat, ani_str, day_idx, s2_w2p_ax, 'usw', 'w2p', 's2_manual_restrict', 0, max_presp, ball_size_range);
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s2_w2p_ax, 'bsw', 'w2p', 's2_manual_restrict', 0, max_presp, ball_size_range);
    %plot_3d_s1s2(all_dat, ani_str, day_idx, s2_w2p_ax, 'mw', 'w2p', 's2_manual_restrict', 1, max_presp, ball_size_range);
    

    %%
    disp('FIX : w1p/w2p labels ; line for x/y length legend ; add dummy balls for legend ');
   % print_fig_LR(fh, 'mw_S1_w1p')

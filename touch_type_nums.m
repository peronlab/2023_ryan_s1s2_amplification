function touch_type_nums 
    ani_list= {'an014351', 'an017518', 'an014359', 'an016623', 'an017517', 'an016650','an016652','an014363','an017510'};
    date_list_1= {'06_26', '08_23', '08_25', '07_29', '08_22', '05_27','06_03','08_08','08_17'};
    date_list_2= {'06_27', '08_24', '08_26', '08_02', '08_23', '05_31','06_08','08_09','08_19'}; 
    root_dir= '/Volumes/volume_imaging/';
    all_dat = get_s1s2_all_dat;
    ani_used = find(all_dat.has_deep_data_s1s2);

    
    for i=1:length(ani_list)
        
        ai= ani_used(i);
        anim=string(ani_list(i));
        cd(sprintf('%s/%s/session_neuropilone/',root_dir, anim)); %go to animal directory 
        load (sprintf('%s_2022_%s_sess.mat', anim, string(date_list_1(i))))
        wdat_day1= s.getWhiskerTouchTrialTypes; 
        %load (sprintf('%s_2022_%s_sess.mat', anim, string(date_list_2(i))))
        %wdat_day2= s.getWhiskerTouchTrialTypes; 
        
        n_type (i, 1)= length(wdat_day1.w1ExclusiveTouchTrials)%+ %length(wdat_day2.w1ExclusiveTouchTrials); %w1ex
        n_type (i, 2)= length(intersect(wdat_day1.w1ProTouchTrials, wdat_day1.w1ExclusiveTouchTrials))%+ length(intersect(wdat_day2.w1ProTouchTrials, wdat_day2.w1ExclusiveTouchTrials)); %w1p 
        n_type (i, 3)= length(intersect(wdat_day1.w1RetTouchTrials, wdat_day1.w1ExclusiveTouchTrials))%+ length(intersect(wdat_day2.w1RetTouchTrials, wdat_day2.w1ExclusiveTouchTrials)); %w1r
        n_type (i, 4)= length(wdat_day1.w2ExclusiveTouchTrials)%+ length(wdat_day2.w2ExclusiveTouchTrials); %w2ex
        n_type (i, 5)= length(intersect(wdat_day1.w2ProTouchTrials, wdat_day1.w2ExclusiveTouchTrials))%+ length(intersect(wdat_day2.w2ProTouchTrials, wdat_day2.w2ExclusiveTouchTrials)); %w2p
        n_type (i, 6)= length(intersect(wdat_day1.w2RetTouchTrials, wdat_day1.w2ExclusiveTouchTrials))%+ length(intersect(wdat_day2.w2RetTouchTrials, wdat_day2.w2ExclusiveTouchTrials)); %w2r
        n_type (i, 7)= length(wdat_day1.multiWhiskerTrials)%+ length(wdat_day2.multiWhiskerTrials); %multi_all

        w2p_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w2ProTouchTrials); %w2p_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w2ProTouchTrials);
        w1p_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w1ProTouchTrials); %w1p_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w1ProTouchTrials);
        w2r_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w2RetTouchTrials); %w2r_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w2RetTouchTrials);
        w1r_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w1RetTouchTrials); %w1r_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w1RetTouchTrials);

        n_type (i, 8)= length(intersect(w2p_d1,w1p_d1))%+ length(intersect(w2p_d2,w1p_d2));%w2p/w1p
        n_type (i, 9)= length(intersect(w2r_d1,w1p_d1))%+ length(intersect(w2r_d2,w1p_d2));%w2r/w1p
        n_type (i, 10)= length(intersect(w2p_d1,w1r_d1))%+ length(intersect(w2p_d2,w1r_d2)); %w2p/w1r
        n_type (i, 11)= length(intersect(w2r_d1,w1r_d1))%+ length(intersect(w2r_d2,w1r_d2)); %w2r/w1r
        
        number_trials (i) = (length(wdat_day1.allTrials)-length(wdat_day1.noTouchTrials))/length(wdat_day1.allTrials); %this is fraction of trials with touch
        
        num_cells(i)= size((all_dat.anim_mats{ai}.probRespW1),2);
    end 
    
        mean_plot= mean(n_type,1);
%% this is for example animal, for now, using all anims. 
%     load an016652_2022_06_01_sess.mat;
%     wdat_day1=  s.getWhiskerTouchTrialTypes;
%     %%
%     load an016652_2022_06_03_sess.mat;
%     wdat_day2=  s.getWhiskerTouchTrialTypes;
    %% order: w1ex., w1p(ex), w1r(ex), w2ex, w2p (ex), w2r(ex),multi(w1w2all), w2pw1p, w2rw1p, w1pw2r, w1r2r
%     n_type (1)= length(wdat_day1.w1ExclusiveTouchTrials)+ length(wdat_day2.w1ExclusiveTouchTrials); %w1ex
%     n_type (2)= length(intersect(wdat_day1.w1ProTouchTrials, wdat_day1.w1ExclusiveTouchTrials))+ length(intersect(wdat_day2.w1ProTouchTrials, wdat_day2.w1ExclusiveTouchTrials)); %w1p 
%     n_type (3)= length(intersect(wdat_day1.w1RetTouchTrials, wdat_day1.w1ExclusiveTouchTrials))+ length(intersect(wdat_day2.w1RetTouchTrials, wdat_day2.w1ExclusiveTouchTrials)); %w1r
%     n_type (4)= length(wdat_day1.w2ExclusiveTouchTrials)+ length(wdat_day2.w2ExclusiveTouchTrials); %w2ex
%     n_type (5)= length(intersect(wdat_day1.w2ProTouchTrials, wdat_day1.w2ExclusiveTouchTrials))+ length(intersect(wdat_day2.w2ProTouchTrials, wdat_day2.w2ExclusiveTouchTrials)); %w2p
%     n_type (6)= length(intersect(wdat_day1.w2RetTouchTrials, wdat_day1.w2ExclusiveTouchTrials))+ length(intersect(wdat_day2.w2RetTouchTrials, wdat_day2.w2ExclusiveTouchTrials)); %w2r
%     n_type (7)= length(wdat_day1.multiWhiskerTrials)+ length(wdat_day2.multiWhiskerTrials); %multi_all
%     
%     w2p_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w2ProTouchTrials); w2p_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w2ProTouchTrials);
%     w1p_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w1ProTouchTrials); w1p_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w1ProTouchTrials);
%     w2r_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w2RetTouchTrials); w2r_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w2RetTouchTrials);
%     w1r_d1= intersect(wdat_day1.multiWhiskerTrials, wdat_day1.w1RetTouchTrials); w1r_d2= intersect(wdat_day2.multiWhiskerTrials, wdat_day2.w1RetTouchTrials);
%     
%     n_type (8)= length(intersect(w2p_d1,w1p_d1))+ length(intersect(w2p_d2,w1p_d2));%w2p/w1p
%     n_type (9)= length(intersect(w2r_d1,w1p_d1))+ length(intersect(w2r_d2,w1p_d2));%w2r/w1p
%     n_type (10)= length(intersect(w2p_d1,w1r_d1))+ length(intersect(w2p_d2,w1r_d2)); %w2p/w1r
%     n_type (11)= length(intersect(w2r_d1,w1r_d1))+ length(intersect(w2r_d2,w1r_d2)); %w2r/w1r
%     
    labels= {'w1. ex all', 'w1p', 'w1r', 'w2ex all', 'w2p', 'w2r', 'multi_all', 'w2p/w1p', 'w2r/w1p', 'w2p/w1r', 'w2r/w1r'}
    %% 
    fh= figure;
    ax=axes; 
    x= [2 3 5 6 8]; 
    sub_sel= [2 3 5 6 7];
    
    bar(ax, x(1:2), mean_plot(2:3),'r'); hold on;
    plot (x, n_type(:,sub_sel),'ok');
    bar (ax, x(3:4), mean_plot(5:6),'b'); 
    bar (ax, x(5), mean_plot(7), 'k');
    %x=1:12;
    y=150;
    plot(x,y*ones(length(x)))
    %ax.YLim=[0 150]
    ax.TickDir='out';ax.Box='off';ax.XTick=[2 3 5 6 8];set(ax,'XTickLabel', labels(sub_sel)); 
    %print_fig_LR(fh, '
    
    %% num touches calc 
    mean_touchtrialfrac= mean(number_trials); 
    SEM= std(number_trials)/sqrt(length(number_trials));
    mean_num_cells= mean(num_cells); 
    SEM_nc= std(num_cells)/sqrt(length(num_cells));
    
end 
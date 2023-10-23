% testing a figure that takes the mean response of each cell, adds it all
% up, and then divides by category (gives as fraction) 

    fh = figure('Position',[0 0 1500 900]);
    ax= axes;

    all_dat = get_s1s2_all_dat;

    %s1_tag = 's1_manual restrict';
    %s1_tag = 's1_manual_restrict';
    %s2_tag = 's2_all';
    %s2_tag = 's2_manual_restrict';
    vali= find(all_dat.has_deep_data_s1s2); 

    %restrict_tag= s1_tag; 
    
    cts = {'usw', 'bsw', 'mw'};

    ratio_mat = zeros(3,length(vali));
 for aa=1:2
     
    for v=1:length(vali)
        ai=vali(v); 
        
     if aa== 1
         % restrict to desired area's id range
        area_nrni = get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, 's1_manual_restrict', all_dat, ai);
        d_s1_used = 1;
        if (strcmp(all_dat.anims{ai}, 'an016652'))
            d_s1_used=3;
        end 
     end 
     if aa==2 
        area_nrni = get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, 's2_manual_restrict', all_dat, ai);
        d_s1_used = 3;
        if (strcmp(all_dat.anims{ai}, 'an016652')) || (strcmp(all_dat.anims{ai}, 'an014359'))
            d_s1_used=1;
        end 
     end 
     
     
        %di = min(find(all_dat.retro_dates{ai}));
        dat_dayi = d_s1_used;

        % and now get the ratio of actual to expected fraction by type
        for c=1:length(cts) 
            ct = cts{c};
            cti = intersect(area_nrni,all_dat.types_by_idx{ai}.([ct '_by_day']){dat_dayi});
            N_ct(c) = length(cti);
            %N_red_ct(c) = length(intersect(cti,all_dat.types_by_idx{ai}.is_red));
        end
        N_touch = sum(N_ct);
        N_all = length(area_nrni);
       % N_red = length(all_dat.types_by_idx{ai}.is_red);

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
%       mu_resp_per_cell = resp_dff_mat;   

        net_resp_overall = nansum(mu_resp_per_cell);
        
        for ctt= 1:length(cts)
            ct = cts{ctt};
            cti = intersect(area_nrni,all_dat.types_by_idx{ai}.([ct '_by_day']){dat_dayi});
            mu_resp_per_ID(ctt)= nansum(mu_resp_per_cell(cti));% gives sum per ID type

        end 
        
        frac_resp (:, ai,aa)= mu_resp_per_ID/net_resp_overall;
    end 
 end 

 %%
    frac_of_resp= frac_resp(:,vali,:); 
    ms = 10;
    
    %hold(ax, 'on');
    %title(ax, tstr);

    %xoffs = rand(size(ratio_mat,2),1)*0.2-0.1;
    %xoffs = xoffs*0; % disable jitter
    
    hold(ax,'on')
    x1= [1.5 2.5 3.5]; %plot
    x2= [4.5 5.5 6.5];
    color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
% 
   for cts=1:3   
       bar (ax, x1(cts),mean((frac_of_resp(cts,:,1)),'omitnan'),'FaceColor',color_arr(cts,:)); hold on; %s1
       plot (ax, x1(cts),((frac_of_resp(cts,:,1))),'ko')
       %title (ax, tstr); 
       set(ax,'TickDir','out','FontSize',15);
       bar (ax, x2(cts),mean((frac_of_resp(cts,:,2)),'omitnan'),'FaceColor',color_arr(cts,:)); hold on; %s1
       plot (ax, x2(cts),((frac_of_resp(cts,:,2))),'ko')
   end 
    

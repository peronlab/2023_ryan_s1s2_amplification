%select an example animal and get correlation matrices

%function corr_mats_example_ani
    
    fig_type=3 ; %change if not for figure 2, will change an nums etc. figure 2 this number should be 2

    all_dat= get_s1s2_all_dat; 
    
    if fig_type==2 
        ani_used= find(all_dat.has_deep_data_s1s2);
        ani= ani_used(2); %usable in this code: 2 4 5 8 9
        day_use=2;
        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{ani}.ids, 's2_manual_restrict', all_dat, ani);
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{ani}.ids, 's1_manual_restrict', all_dat, ani);

    else 
        ani_used= find(all_dat.retrograde_label_type==1);
        ani= ani_used(3); 
        day_use=1;
        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{ani}.ids, 's2_manual_restrict', all_dat, ani);
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{ani}.ids, 's1_manual_restrict', all_dat, ani);
        s2i= intersect (s2i, all_dat.types_by_idx{ani}.is_red);
        s1i= intersect (s1i, all_dat.types_by_idx{ani}.is_red); %for fig 3 we wanna restrict to proj. neurons 
    end 
    
    %% 
    date_str= all_dat.valid_dates{ani}{day_use}; %Find the correct date to use as a string for file opening 
    anim= all_dat.anims{ani}; %also find the correct anim string for file opening
    
    
    
        
     
    cd(sprintf('%s/%s/session_neuropilone/pairwiseCorrelations/',all_dat.root_dir, anim)) %go to animal directory 
    
    if fig_type==2 
        sv_str='001'; %str name for the subvolume using (this is what the most superficial s1 is usually called
        sv_str_s2= '091';

        s1= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str)); %dat, load the correct file as DAT 
        s2= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str_s2)); 
    else 
        %sv_str= '001'; 
        sv_str= '091';
        %s1= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str));
        %s2= load (sprintf('%s_2023_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str));

    end 
        cts = {'usw', 'bsw', 'mw'};
        non_merged = find(~all_dat.anim_data{ani}.is_merged_across_days);
        day_used_S1= non_merged(day_use);
    
    for c=1:3 
        a=ani; 
        ct= cts{c};
        if fig_type==2
        cti = intersect(s1i,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1});
        %idi = find(ismember(s1.dat.cellIds, all_dat.anim_data{a}.ids(cti)));
        %cti_corrs= s1.dat.corrMat(idi,idi); %gives corr matrix for touch type
        %cti_corrs_reflect= triu(cti_corrs)+(triu(cti_corrs,1)'); %mirrors corr mat about the diagonal to more easily compute
        
        %cti_corrs=sort(cti_corrs,2,'descend');
                for tt=1:1 % loop through 4 touch types: w1p, w1r, w2p, w2r. trying to find s1, mw cells that specifically respond to each touch type 
                    if tt==1 touch_str= 'whiskerW1PExcTouchTrials'; type='ev_touch_w1p_by_day'; end 
                    if tt==2 touch_str= 'whiskerW2PExcTouchTrials'; type='ev_touch_w2p_by_day'; end 
                    if tt==3 touch_str= 'whiskerW1RExcTouchTrials'; type='ev_touch_w1r_by_day'; end 
                    if tt==4 touch_str= 'whiskerW2RExcTouchTrials'; type='ev_touch_w2r_by_day';end 
               
                    yr_str=2022; 
                    ia= load (sprintf('%s_2022_%s_sess__sv_%s_%s.mat',anim,  date_str, sv_str, touch_str));
                    resp_cells=intersect(cti, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, s1i)));
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    %num_cell (c, tt, ai)= length(idi); %sanity check if we need cell # cutoffs
                    cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type --> these are correlations with OTHER cells of the exact same type
                    cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                    %cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                    
                    %clear idi; clear resp_cells; 
                    
                end 
        
       %%
                
        %%
        cti_s2= (intersect(s2i,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}));%-(min(s2i-1));
        %idi2 = find(ismember(s2.dat.cellIds, all_dat.anim_data{a}.ids(cti_s2)));
        %cti_corrs_s2= (s2.dat.corrMat(idi2,idi2));
        %cti_corrs_reflect_s2= (triu(cti_corrs_s2))+(triu(cti_corrs_s2,1)');
        
              for tt=1:1 % loop through 4 touch types: w1p, w1r, w2p, w2r. trying to find s1, mw cells that specifically respond to each touch type 
                    if tt==1 touch_str= 'whiskerW1PExcTouchTrials'; type='ev_touch_w1p_by_day'; end 
                    if tt==2 touch_str= 'whiskerW2PExcTouchTrials'; type='ev_touch_w2p_by_day'; end 
                    if tt==3 touch_str= 'whiskerW1RExcTouchTrials'; type='ev_touch_w1r_by_day'; end 
                    if tt==4 touch_str= 'whiskerW2RExcTouchTrials'; type='ev_touch_w2r_by_day';end 
               
                    
                    ia2= load (sprintf('%s_2023_%s_sess__sv_%s_%s.mat',anim, date_str, sv_str_s2, touch_str));
                    resp_cells2=intersect(cti_s2, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, s2i)));
                    idi2= find(ismember(ia2.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells2))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    %num_cell (c, tt, ai)= length(idi); %sanity check if we need cell # cutoffs
                    cti_corrs_s2= ia2.dat.corrMat(idi2,idi2); %gives corr matrix for touch type --> these are correlations with OTHER cells of the exact same type
                    cti_corrs_reflect_s2= triu(cti_corrs_s2)+triu(cti_corrs_s2,1)'; %mirrors corr mat about the diagonal to more easily compute
                    %cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                    
                end 
        
        
        %sort s1 by max in each col
        [t,idx]= intersect(max(cti_corrs_reflect), sort(max(cti_corrs_reflect),'descend'));
        idx=flip(idx); 
        cti_corrs_reflect= cti_corrs_reflect(idx,idx);
        
        %%% sort s2 by max in each col
        [t,idx2]= intersect(max(cti_corrs_reflect_s2), sort(max(cti_corrs_reflect_s2),'descend'));
        idx2=flip(idx2); 
        cti_corrs_reflect_s2= cti_corrs_reflect_s2(idx2,idx2);
        
        
        if c==1 usw_S1= cti_corrs_reflect; usw_S2=cti_corrs_reflect_s2; end 
        if c==2 bsw_S1= cti_corrs_reflect; bsw_S2=cti_corrs_reflect_s2; end 
        if c==3 mw_S1= cti_corrs_reflect; mw_S2=cti_corrs_reflect_s2;end 
        
        clear cti_corrs_reflect; clear cti_corrs_reflect_s2; 
        
        else
        %s2i= s1i;
        cti = intersect(s2i,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1});
        %idi = find(ismember(s1.dat.cellIds, all_dat.anim_data{a}.ids(cti)));
        %cti_corrs= s1.dat.corrMat(idi,idi); %gives corr matrix for touch type
        %cti_corrs_reflect= triu(cti_corrs)+(triu(cti_corrs,1)'); %mirrors corr mat about the diagonal to more easily compute

        for tt=3:3% loop through 4 touch types: w1p, w1r, w2p, w2r. trying to find s1, mw cells that specifically respond to each touch type 
                    if tt==1 touch_str= 'whiskerW1PExcTouchTrials'; type='ev_touch_w1p_by_day'; end 
                    if tt==2 touch_str= 'whiskerW2PExcTouchTrials'; type='ev_touch_w2p_by_day'; end 
                    if tt==3 touch_str= 'whiskerW1RExcTouchTrials'; type='ev_touch_w1r_by_day'; end 
                    if tt==4 touch_str= 'whiskerW2RExcTouchTrials'; type='ev_touch_w2r_by_day';end 
               
                    %yr_str=2022; 
                    ia= load (sprintf('%s_2023_%s_sess__sv_%s_%s.mat',anim,  date_str, sv_str, touch_str));
                    resp_cells=intersect(cti, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, s2i)));
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    %num_cell (c, tt, ai)= length(idi); %sanity check if we need cell # cutoffs
                    cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type --> these are correlations with OTHER cells of the exact same type
                    cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                    %cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                    
                    %clear idi; clear resp_cells; 
                    
                end 
        
        %sort s1 by max in each col
        [t,idx]= intersect(max(cti_corrs_reflect), sort(max(cti_corrs_reflect),'descend'));
        idx=flip(idx); 
        cti_corrs_reflect= cti_corrs_reflect(idx,idx);
        
        
        
        if c==1 usw_S1= cti_corrs_reflect; end %usw_S2=cti_corrs_reflect_s2; end 
        if c==2 bsw_S1= cti_corrs_reflect; end %bsw_S2=cti_corrs_reflect_s2; end 
        if c==3 mw_S1= cti_corrs_reflect;  end %mw_S2=cti_corrs_reflect_s2;end  
        
       
            
        end 
    end 
    
    if fig_type==2 
    maxi=0.3;
    cm=colormap_human(64);
    figure; heatmap(usw_S1, 'Colormap', cm,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'USW S1'); 
    figure; heatmap(usw_S2, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'USW S2');
    figure; heatmap(bsw_S1, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'bSW S1');
    figure; heatmap(bsw_S2, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'bSW S2');
    figure; heatmap(mw_S1, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'mw S1');
    figure; heatmap(mw_S2, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'mw S2');
    else 
    maxi=0.3;
    cm=colormap_human(64);
    figure; heatmap(usw_S1, 'Colormap', cm,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'USW S1'); 
    %figure; heatmap(usw_S2, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'USW S2');
    figure; heatmap(bsw_S1, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'bSW S1');
    %figure; heatmap(bsw_S2, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'bSW S2');
    figure; heatmap(mw_S1, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'mw S1');
    %figure; heatmap(mw_S2, 'Colormap', jet,'GridVisible', 'off', 'ColorLimits', [0 maxi], 'Title', 'mw S2');
    end 
    
%end 


%% test

%  diag_sort= sort(diag_val, 'descend');
%         [share_val, idx]= intersect(diag_val,diag_sort);
%         idx=flip(idx); 
%         index_use= idi(idx');
%         cti_corrs= s1.dat.corrMat(index_use,index_use);
%         
%         S1_corrs=s1.dat.corrMat; S1_corrs= triu(S1_corrs)+(triu(S1_corrs,1)');
%         S1_sort=S1_corrs(index_use, index_use);
%         
%         for i=1:(length(cti_corrs)-1)
%              diag_val(i)= cti_corrs(i,i+1);    
%         end

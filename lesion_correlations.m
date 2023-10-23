%% function for finding correlations pre and post lesion 

%%
%function find_lesion_corrs
all_dat=get_s1s2_all_dat;
s1_lesion_ani = find(~strcmp('',all_dat.s1_lesion_date)); 
%s1_lesion_ani= s1_lesion_ani([1:4, 6:end]);%20048 has degenerate spontaneous correlations, removed for now
s2_lesion_ani = find(~strcmp('',all_dat.s2_lesion_date));
sham_lesion_ani = find(~strcmp('',all_dat.sham_lesion_date)); % use all shams

[pre_S2_les, post_S2_les] = lesion_corrs(all_dat, s2_lesion_ani, all_dat.lesion_dates, 'S2 lesion',1) %0: stim locked, 1: spont,2: all cells, not w.in group 

%end 
%% 
function [stim_corr_type_pre, stim_corr_type_post]= lesion_corrs(all_dat,ani_used, dates_used, lesion_type, spont)   
%% setup 
    fh = figure('Position', [0 0 800 800]); 

    ax_pre= subplot('Position', [.1 .1 0.3 0.8]);  
    hold on;
    ax_post= subplot('Position', [.5 .1 0.3 0.8]);   
    hold on;

    for ai=1:length(ani_used)
        a = ani_used(ai);
        
        pre_di= find(dates_used{a}==-1); %find pre lesion day. 
        pre_di= pre_di(1); %if there is more than one day choose the last. 
        post_di= find(dates_used{a}==1); %find post lesion day 
        post_di= post_di(1); %if there is more than one day choose the first. 
        
        if lesion_type=='S2 lesion'
           sv_str= 'none'%'001'; %if you lesion s2, you are imaging in s1 pre/post. most lesion guys dont have an sv str for the pre/post day
           area_ofi= 'S1';
           si = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);
           
        end 
        
        if lesion_type== 'S1 lesion'
            sv_str= '091'; %str name for the subvolume using 
            area_ofi= 'S2'; 
            si = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
        end
        
        if (strcmp(all_dat.anims{a}, 'an016623'))
            sv_str='001';    
        end 
        
        if (strcmp(all_dat.anims{a}, 'an014351'))
            sv_str='101';
        end 
         
        if(strcmp(all_dat.anims{a}, 'an018927')) || (strcmp(all_dat.anims{a}, 'an020048')) || (strcmp(all_dat.anims{a}, 'an018922'))
            sv_str= '001';
        end 

        yr_str= '2022';% set year to load to 2023, but change to 2022 for specific animals below
        
        anim= all_dat.anims{a}; %also find the correct anim string for file opening
        
        cd(sprintf('%s/%s/session_neuropilone/pairwiseCorrelations/',all_dat.root_dir, anim)) %go to animal directory 
       
        
        %ia= load (sprintf('%s_%s_%s_sess__sv_%s_allTime.mat',anim, yr_str, date_str, sv_str)); %dat, load the correct file as DAT 
   
        %start by finding MW Cells in one animal. day 1
        non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);
        day_used_S1_pre= non_merged(pre_di);
        cts = {'usw', 'bsw', 'mw'};
             
            % NOTE: right now were in one subvolume only, will probably
            % need both. 
            %% find correlations for each touch cell type and each single whisker touch type
       if spont==0 
       for dd=1:2 %go through pre and post dates
           if dd==1 
               day_used_S1= non_merged(pre_di); 
               date_str=all_dat.valid_dates{a}{pre_di}; 
               if (strcmp(all_dat.anims{a}, 'an018933'))
                sv_str='001'; %mistake in labeling this one guy's subvolume
               end 
           end 
           if dd==2 day_used_S1= non_merged(post_di); date_str= all_dat.valid_dates{a}{post_di};
                if (strcmp(all_dat.anims{a}, 'an018933'))
                sv_str='none'; %mistake in labeling this one guy's subvolume
                end 
               if (strcmp(all_dat.anims{a}, 'an018927')) || (strcmp(all_dat.anims{a}, 'an020048'))||(strcmp(all_dat.anims{a}, 'an018922'))
                   sv_str= 'none'; 
               end 
           end 
           
             for c=1:length(cts) 
                ct = cts{c};
                cti = intersect(si,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}); %gives cells in area ID range that are mW, now need to break up further by touch type
              

                for tt=1:4 % loop through 4 touch types: w1p, w1r, w2p, w2r. trying to find s1, mw cells that specifically respond to each touch type 
                    if tt==1 touch_str= 'whiskerW1PExcTouchTrials'; type='ev_touch_w1p_by_day'; end 
                    if tt==2 touch_str= 'whiskerW2PExcTouchTrials'; type='ev_touch_w2p_by_day'; end 
                    if tt==3 touch_str= 'whiskerW1RExcTouchTrials'; type='ev_touch_w1r_by_day'; end 
                    if tt==4 touch_str= 'whiskerW2RExcTouchTrials'; type='ev_touch_w2r_by_day';end 
               
                    if strcmp(sv_str,'none')
                        ia= load (sprintf('%s_%s_%s_sess__%s.mat',anim, yr_str, date_str, touch_str));
                    else  
                        ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, sv_str, touch_str));
                    end 
                    resp_cells=intersect(cti, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, si)));
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type --> these are correlations with OTHER cells of the exact same type
                    cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                    cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                    mean_corr_AT(c,tt)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type and for each touch type in area of interest
                end 
               
             end   
         
             
            means= nanmean(mean_corr_AT,2);
            mean_plot_rlly(:,ai,dd)= means; %columns are animals, rows are cell type: usw/bsw/mw. and dd 1 or 2 is pre/post
            
            
       end 
       end
    
       if spont==1 %if youre using spontaneous instead of stim locked, last part of function 
            for dd=1:2 %go through pre and post dates
           if dd==1 
               day_used_S1= non_merged(pre_di); 
               date_str=all_dat.valid_dates{a}{pre_di}; 
               if (strcmp(all_dat.anims{a}, 'an018933'))
                sv_str='001'; %mistake in labeling this one guy's subvolume
               end 
           end 
           if dd==2 day_used_S1= non_merged(post_di); date_str= all_dat.valid_dates{a}{post_di};
                if (strcmp(all_dat.anims{a}, 'an018933'))
                sv_str='none'; %mistake in labeling this one guy's subvolume
                end 
               if (strcmp(all_dat.anims{a}, 'an018927')) || (strcmp(all_dat.anims{a}, 'an020048'))||(strcmp(all_dat.anims{a}, 'an018922'))
                   sv_str= 'none'; 
               end 
           end 
           
             for c=1:length(cts) 
                ct = cts{c};
                cti = intersect(si,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}); %gives cells in area ID range that are mW, now need to break up further by touch type
                touch_str= 'whiskerNoTouchTrials';
               
                    if strcmp(sv_str,'none')
                        ia= load (sprintf('%s_%s_%s_sess__%s.mat',anim, yr_str, date_str, touch_str));
                    else  
                        ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, sv_str, touch_str));
                    end 

                    
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(cti)));
                    cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type
                    cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                    cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                    mean_corr_AT(c)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type, s1
                    %mean_corr_ancurr (c,dd)= mean(cti_corrs_per_cell);
                   
                    %mean_corr_AT(c,tt)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type and for each touch type in area of interest
                 
             end   
         
             
            means= (mean_corr_AT);
            mean_plot_rlly(:,ai,dd)= means; %columns are animals, rows are cell type: usw/bsw/mw. and dd 1 or 2 is pre/post
            
            
       end 
       end 
       
       if spont==2 
           for dd=1:2 
               if dd==1 
               day_used_S1= non_merged(pre_di); 
               date_str=all_dat.valid_dates{a}{pre_di}; 
               if (strcmp(all_dat.anims{a}, 'an018933'))
               sv_str='001'; %mistake in labeling this one guy's subvolume
               end 
           end 
           if dd==2 day_used_S1= non_merged(post_di); date_str= all_dat.valid_dates{a}{post_di};
                if (strcmp(all_dat.anims{a}, 'an018933'))
                sv_str='none'; %mistake in labeling this one guy's subvolume
               end 
           end 
            
            
         for c=1:length(cts) 
                ct = cts{c};
                cti = intersect(si,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}); %gives cells in area ID range that are mW, now need to break up further by touch type
              

                for tt=1:4 % loop through 4 touch types: w1p, w1r, w2p, w2r. trying to find s1, mw cells that specifically respond to each touch type 
                    if tt==1 touch_str= 'whiskerW1PExcTouchTrials'; type='ev_touch_w1p_by_day'; end 
                    if tt==2 touch_str= 'whiskerW2PExcTouchTrials'; type='ev_touch_w2p_by_day'; end 
                    if tt==3 touch_str= 'whiskerW1RExcTouchTrials'; type='ev_touch_w1r_by_day'; end 
                    if tt==4 touch_str= 'whiskerW2RExcTouchTrials'; type='ev_touch_w2r_by_day';end 
               
                    if strcmp(sv_str,'none')
                        ia= load (sprintf('%s_%s_%s_sess__%s.mat',anim, yr_str, date_str, touch_str));
                    else  
                        ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, sv_str, touch_str));
                    end 
                    resp_cells=intersect(cti, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, si)));
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    all_resp_cells= (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, si));
                    all_cells_type_idi= find(ismember(ia.dat.cellIds,all_dat.anim_data{a}.ids(all_resp_cells)));
                    cti_corrs= ia.dat.corrMat(all_cells_type_idi,all_cells_type_idi); %gives corr matrix for touch type --> these are correlations with OTHER cells of the exact same type
                    cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                    idx_to_parse= find (ismember(all_cells_type_idi,idi));
                    cti_corrs_per_cell= nanmean(cti_corrs_reflect(:,idx_to_parse),1); %gives mean correlation of each cell with every other cell that responds to same type
                    mean_corr_AT(c,tt)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type and for each touch type in area of interest
                   
                end 
         end
         
            means= nanmean(mean_corr_AT,2);
            mean_plot_rlly(:,ai,dd)= means; %columns are animals, rows are cell type: usw/bsw/mw. and dd 1 or 2 is pre/post
            
       end 
    end 
    end 
     x1= [1.5 3.5 5.5 7.5]; %x2= [2.5 4.5 6.5 8.5];   
     color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color ; 0.75 0.75 0.75];
     names={'usw', 'bsw', 'mw'};
    
   
     for d=1:2
         if d==1 ax=ax_pre; end 
         if d==2 ax=ax_post; end 
         
         for ctns=1:3 
             n(ctns)= sum(~isnan(mean_plot_rlly(ctns,:,d)));
             bar (ax, x1(ctns),nanmean(mean_plot_rlly(ctns,:,d)),'FaceColor',color_arr(ctns,:)); hold on; %s1
             plot (ax, x1(ctns),(nonzeros(mean_plot_rlly(ctns,:,d))),'ko','MarkerSize',10)
             ax.TickDir='out';ax.Box='off';
             ax.YLim= [0 0.15]
         end  
     end 
      
       string_title= sprintf('stim locked corrs pre/post lesion. imaging in: %s, usw/bsw/mw n= %i/%i/%i', area_ofi, n(1),n(2),n(3))
       title (ax, string_title)
        
       stim_corr_type_pre=mean_plot_rlly(:,:,1);
       stim_corr_type_post= mean_plot_rlly(:,:,2);
     
end 

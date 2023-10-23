%% function for spontaneous correlations in projecting neurons 

function [stim_corr_type]= spont_corrs_better_sub(retro, proj_only,sub_better,discard_an, ylimit)   
%% setup 

    all_dat = get_s1s2_all_dat;

    fh = figure('Position', [0 0 500 500]); 

    ax= subplot('Position', [.1 .1 0.8 0.8]);   
    hold on;
    
    ani_use= find(all_dat.retrograde_label_type==retro);
    ani_used=setdiff(ani_use, discard_an) %take out animals with not enough cells for one of the categories
    
    idx=find(ismember(ani_use, discard_an));
    array_ids=1:length(ani_use);
    id_use= setdiff(array_ids,idx);
    sub_better= sub_better(id_use)
    
    %% pulling correlations
    %% 
    for ai=1:length(ani_used)
        a = ani_used(ai);
        day_use= 1; %for now, just pick the 2nd pre lesion day for all these guys, and use the superficial subvolume. 
        
        if retro==1
           sv_str(1,:)= '091';
           %sv_str= '091';
           sv_str(2,:)= '092';
           id_max= 920000000;
           area_ofi= 'S2'
           si = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
           
        end 
        
        if retro==2 
            sv_str (1,:)='001'; %str name for the subvolume using (this is what the most superficial s1 is usually called
            sv_str (2,:)='002';
            id_max= 20000000;
            area_ofi= 'S1'; 
            si = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);
        end
        
        
        
        yr_str= '2023';% set year to load to 2023, but change to 2022 for specific animals below
        
        
        if (strcmp(all_dat.anims{a}, 'an014359')) ||(strcmp(all_dat.anims{a}, 'an014355')) ||(strcmp(all_dat.anims{a}, 'an014354')) ||(strcmp(all_dat.anims{a}, 'an014362')) ||(strcmp(all_dat.anims{a}, 'an018912')) ||(strcmp(all_dat.anims{a}, 'an018919')) || (strcmp(all_dat.anims{a}, 'an020055'))  
            yr_str= '2022'; 
        end 
        
        if (strcmp(all_dat.anims{a}, 'an020035'))
            sv_str(1,:)='001'; %mistake in labeling this one guy's subvolume
            sv_str (2,:)='002';
        end 
        
        
        date_str= all_dat.valid_dates{a}{day_use}; %Find the correct date to use as a string for file opening 
        anim= all_dat.anims{a}; %also find the correct anim string for file opening
        
         if proj_only==1
             si= intersect(si,all_dat.types_by_idx{a}.is_red);
         end 

         if proj_only==0
             si= setdiff(si,(intersect(si,all_dat.types_by_idx{a}.is_red)));
         end   
         
         if sub_better(ai)==3
             subvol_1= sv_str(1,:); 
             subvol_2= sv_str(2,:);
         else 
            subvol= sv_str(sub_better(ai),:);

         end 
         
         if (strcmp(all_dat.anims{a}, 'an014359')) || (strcmp(all_dat.anims{a}, 'an014362'))|| (strcmp(all_dat.anims{a}, 'an014355'))|| (strcmp(all_dat.anims{a}, 'an014354')) %old animals with bigger dzs, just use one subvol
                subvol= sv_str(1,:);
                subvol_1= sv_str(1,:);
                subvol_2= sv_str(1,:);
         end 
         
         
         
        cd(sprintf('%s/%s/session_neuropilone/pairwiseCorrelations/',all_dat.root_dir, anim)) %go to animal directory 
        touch_str= 'allTime';
       
        %ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, subvol, touch_str));

        %ia= load (sprintf('%s_%s_%s_sess__sv_%s_allTime.mat',anim, yr_str, date_str, sv_str)); %dat, load the correct file as DAT 
   
        %start by finding MW Cells in one animal. day 1
         non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);
         day_used_S1= non_merged(day_use);
         cts = {'usw', 'bsw', 'mw'};
             
            % NOTE: right now were in one subvolume only, will probably
            % need both. 
            %% find correlations for each touch cell type and each single whisker touch type
      
         if sub_better(ai)==3 
            for dd=1:2 
                if dd==1
                    subvol= subvol_1;
                    ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, subvol, touch_str));

                else 
                    subvol= subvol_2; 
                    ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, subvol, touch_str));

                end 
             for c=1:length(cts) 
            
                ct = cts{c};
                cti = intersect(si,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}); %gives cells in area ID range that are mW, now need to break up further by touch type
                idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(cti)));
                cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type
                cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                mean_corr_AT(c,dd,ai)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type, s1
                mean_corr_ancurr (c,dd)= mean(cti_corrs_per_cell);
                 
                
                
             end   
         
            end 
            means= cat(2,mean_corr_ancurr(:,1),mean_corr_ancurr(:,1));
            mean_plot_rlly(:,ai)= nanmean(means,2); 
            end 
            
           
            
         
         
         if sub_better(ai)==1 || sub_better(ai)==2
            ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, subvol, touch_str));

            for c=1:length(cts) 
                ct = cts{c};
                cti = intersect(si,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}); %gives cells in area ID range that are mW, now need to break up further by touch type
                idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(cti)));
                cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type
                cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                mean_corr_AT(c,ai)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type, s1
                mean_corr_ancurr2(c)=mean(cti_corrs_per_cell);
            end 
         
         %means= mean_corr_ancurr(:,:,dd); %means for this animal, for this subvolume!
         mean_plot_rlly (:,ai)= (mean_corr_ancurr2);
         %mean_plot(:,dd,ai)= nanmean(means,2); %columns are animals, rows are cell type: usw/bsw/mw   
    
            %means= cat(2,mean_corr_ancurr(:,:,1),mean_corr_ancurr(:,:,2));

            %mean_plot(:,dd,ai)= nanmean(means,2); %columns are animals, rows are cell type: usw/bsw/mw
            %mean_plot_rlly(:,ai)= nanmean(mean_corr_AT,2);
            %leaving out NT cells rn
            %s_ci = find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(si)));
            %mean_corr_AT(ai,4,1) = nanmean(reshape(ia.dat.corrMat(s_ci,s_ci),[],1)) ; % ALL cells
         end  
    end 
    
      x1= [1.5 3.5 5.5 7.5]; %x2= [2.5 4.5 6.5 8.5];   
     color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color ; 0.75 0.75 0.75];
     names={'usw', 'bsw', 'mw'};
        for ctns=1:3 
         
         %string_title= sprintf('imaging in: %s, Retro: %i, Projection: %i', area_ofi, retro,proj_only)
         %title (ax, string_title)
         bar (ax, x1(ctns),nanmean(nonzeros(mean_plot_rlly(ctns,:))),'FaceColor',color_arr(ctns,:)); hold on; %s1
         plot (ax, x1(ctns),(nonzeros(mean_plot_rlly(ctns,:))),'ko','MarkerSize',10)
         n(ctns)= length((nonzeros(mean_plot_rlly(ctns,:))));
         %bar (ax, x2(cti),nanmean(nonzeros(mean_corr_AT(:,cti,2))),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
         %plot (ax, x2(cti),(nonzeros(mean_corr_AT(:,cti,2))),'ko')
         ax.TickDir='out';ax.Box='off';
         ylim([0 ylimit])

     end  
       string_title= sprintf('spont corrs. imaging in: %s, Retro: %i, Projection: %i, usw/bsw/mw %i/%i/%i', area_ofi, retro, proj_only, n(1), n(2), n(3))
       title (ax, string_title)
        
       stim_corr_type=mean_plot_rlly;
     
end 
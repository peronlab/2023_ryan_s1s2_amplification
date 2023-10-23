function [sub_better, number_cells_sub_type, sum_cells_subvol, discard_these_anims]= find_usable_subvol(retro, proj_only)
    
    %retro=2; proj_only=1; %tmp
    
    all_dat = get_s1s2_all_dat;
    
    ani_used= find(all_dat.retrograde_label_type==retro);
    
    discard_these_anims=[];

    subvol_bad= zeros(2,length(ani_used));
    
    for ai=1:length(ani_used)
        a = ani_used(ai);
        day_use= 1; %for now, just pick the 2nd pre lesion day for all these guys, and use the superficial subvolume. 
        
        if retro==1
           sv_str= '091';
           sv_str_deep= '092';
           area_ofi= 'S2'
           si = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
           
        end 
        
        if retro==2 
            sv_str='001'; %str name for the subvolume using (this is what the most superficial s1 is usually called
            sv_str_deep='002';
            area_ofi= 'S1'; 
            si = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);
        end
        
        yr_str= '2023';% set year to load to 2023, but change to 2022 for specific animals below
        
        
        if (strcmp(all_dat.anims{a}, 'an014359')) ||(strcmp(all_dat.anims{a}, 'an014355')) ||(strcmp(all_dat.anims{a}, 'an014354')) ||(strcmp(all_dat.anims{a}, 'an014362')) ||(strcmp(all_dat.anims{a}, 'an018912')) ||(strcmp(all_dat.anims{a}, 'an018919')) || (strcmp(all_dat.anims{a}, 'an020055'))  
            yr_str= '2022'; 
        end 
        
        if (strcmp(all_dat.anims{a}, 'an020035'))
            sv_str='001'; %mistake in labeling this one guy's subvolume
            sv_str_deep='002';
        end 
        
        date_str= all_dat.valid_dates{a}{day_use}; %Find the correct date to use as a string for file opening 
        anim= all_dat.anims{a}; %also find the correct anim string for file opening
        
         if proj_only==1
             si= intersect(si,all_dat.types_by_idx{a}.is_red);
         end 

         if proj_only==0
             si= setdiff(si,(intersect(si,all_dat.types_by_idx{a}.is_red)));
         end   
      
        cd(sprintf('%s/%s/session_neuropilone/pairwiseCorrelations/',all_dat.root_dir, anim)) %go to animal directory 
       
        
        %ia= load (sprintf('%s_%s_%s_sess__sv_%s_allTime.mat',anim, yr_str, date_str, sv_str)); %dat, load the correct file as DAT 
   
        %start by finding MW Cells in one animal. day 1
         non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);
         day_used_S1= non_merged(day_use);
         cts = {'usw', 'bsw', 'mw'};
             
            % NOTE: right now were in one subvolume only, will probably
            % need both. 
            %% find correlations for each touch cell type and each single whisker touch type
    for dd=1:2
        if dd==1 sv_str=sv_str; end 
        if dd==2 
           if (strcmp(all_dat.anims{a}, 'an014359')) || (strcmp(all_dat.anims{a}, 'an014362'))|| (strcmp(all_dat.anims{a}, 'an014355'))|| (strcmp(all_dat.anims{a}, 'an014354')) %old animals with bigger dzs, just use one subvol
                sv_str= sv_str; 
           else 
               sv_str=sv_str_deep;
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
               
                    
                    ia= load (sprintf('%s_%s_%s_sess__sv_%s_%s.mat',anim, yr_str, date_str, sv_str, touch_str));
                    resp_cells=intersect(cti, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, si)));
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    num_cell (c, tt, dd, ai)= length(idi); %sanity check if we need cell # cutoffs
                    
                end         
         end
         sum_me= num_cell(:,:,dd, ai);
         number_cells_sub_type(:,dd,ai)=sum(sum_me,2); 
    end
        for i=1:3
         if number_cells_sub_type(i,1,ai)> number_cells_sub_type(i,2,ai) 
             greater_mat (i,1)= 1;
             greater_mat(i,2)=0 ;
         else 
             greater_mat (i,1)= 0;
             greater_mat(i,2)=1 ;
         end 
        end
        
     if sum(greater_mat(:,1)) > sum(greater_mat(:,2))
        sub_better(ai)= 1;%col 1 shows which subvolume wins more categories of cells
     else 
        sub_better(ai)= 2;
     end 
     
     
     sum_cells_subvol(:,ai)= sum(number_cells_sub_type(:,:,ai),1);
    
     for j=1:2 %subvolume 1 to 2
         test= number_cells_sub_type([1 3],j,ai);
         if any(test < 5)
             subvol_bad(j,ai)= 1;
         end 
     end 
     
    
     
     if number_cells_sub_type(3,sub_better(ai),ai) < 5|| number_cells_sub_type(1,sub_better(ai),ai) < 5 || sub_better(ai)==1 && subvol_bad(1,ai)==1 || sub_better(ai)==2 && subvol_bad(2,ai)==1 
         discard_these_anims(end+1)=a;
     end
    
    if (subvol_bad(1,ai)+subvol_bad(2,ai))== 0 %if both subvolumes r up to snuff, the sub_better will print out 3 for that animal, which will tell us to average across subvolumes. if only one is selected, use it, otherwise if there are none that pass anim is discarded
        sub_better(ai)= 3; 
    end 

end 
         
end
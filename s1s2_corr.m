 %% FIRST PASS AT CODE FOR CORRELATiONS LR 2/14/23
 % Right now: this code plots the correlations between touch groups in s1
 % and s2 in the 9 animals we have 2 subvolumes for. using only the most
 % superficial sub per area. YOU can toggle between stim_locked_corrs and
 % non touch spont using the stim-locked= 0 or 1
 

function s1s2_corr
    % setup 
    all_dat = get_s1s2_all_dat;
    projection_restricted = 0; %not totally usable, differenet function for this rn. 
    stim_locked=1; %toggle me to do stim locked versus spont. 
   
    ani_used = find(all_dat.has_deep_data_s1s2) ; %for NOW, know all these corr are right/exist
    %ani_used = find(all_dat.retrograde_label_type == 0); %all non retro animals
    %ani_used= [1 3 12 14 20 21 23 24]; 
    if (projection_restricted)
        ani_used = find(all_dat.retrograde_label_type == 1);
        disp('RESTRICTING TO PROJECTION ; disable this to not do that');
        
        %pause;
    end
   
    fh = figure('Position', [0 0 500 500]); 
    
    ax= subplot('Position', [.1 .1 0.8 0.8]);   
    hold on;
    %frac_all_ax = subplot('Position', [.65 .05 .15 .8]);  
    %hold on;
   
    %% this will all have to be done per animal eventually
    %% 
    for ai=1: length(ani_used)
        a = ani_used(ai);
        day_use= 2; %for now, just pick the 2nd pre lesion day for all these guys, and use the superficial subvolume. 
        sv_str='001'; %str name for the subvolume using (this is what the most superficial s1 is usually called
        id_max= 20000000;
         
        %s2 
        sv_str_s2= '091';
        id_max_s2= 920000000;
        %animal naming convention exceptions
        if (strcmp(all_dat.anims{a}, 'an014359')) %this animal is weird, need to use different day
            day_use = 1;    
        end 
        if (strcmp(all_dat.anims{a}, 'an016652')) 
            sv_str= '091';  
            id_max= 920000000;
            sv_str_s2= '901';
            id_max_s2= 9020000000;
        end 
        if (strcmp(all_dat.anims{a}, 'an014351')) ||(strcmp(all_dat.anims{a}, 'an016650')) 
            sv_str_s2= '101';  
            id_max_s2= 1020000000;
        end 

        if (projection_restricted) ||(strcmp(all_dat.anims{a}, 'an020048'))|| (strcmp(all_dat.anims{a}, 'an018922'))||(strcmp(all_dat.anims{a}, 'an014335'))   % in all these animals we use day 1
            day_use=1;
        end

       date_str= all_dat.valid_dates{a}{day_use}; %Find the correct date to use as a string for file opening 
       anim= all_dat.anims{a}; %also find the correct anim string for file opening
         
         %find IDS : because we want the ids ONLY in one subvolume
         %(superficial), and not as dependent onf area, we need to compute id
         %range a bit differently 
         
 %        s1i = find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < id_max);
         %find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s1_id_range{a}(2));
         %s2i
%         s2i= find(all_dat.anim_data{a}.ids > all_dat.s2_id_range{a}(1) & all_dat.anim_data{a}.ids < id_max_s2);

        s2i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's2_manual_restrict', all_dat, a);
        s1i = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, 's1_manual_restrict', all_dat, a);


%          if (projection_restricted)
%              s1i = intersect(s1i,all_dat.types_by_idx{a}.is_red);
%              s2i = intersect(s2i,all_dat.types_by_idx{a}.is_red);
%          end
         
         s1i_deepinc= find(all_dat.anim_data{a}.ids > all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids < all_dat.s1_id_range{a}(2)); %will need to be able to index by this later
         
         cd(sprintf('%s/%s/session_neuropilone/pairwiseCorrelations/',all_dat.root_dir, anim)) %go to animal directory 
       
         if (stim_locked==0)
         s1= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str)); %dat, load the correct file as DAT 
         s2= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str_s2));
         %s1= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str)); %dat, load the correct file as DAT 
         %s2= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str_s2));
         end 
         
%          if (projection_restricted) %seperating so you can just call the appropriately imaged subvolume, only one area in these guys
%             s1= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str)); %dat, load the correct file as DAT 
%             s2= load (sprintf('%s_2022_%s_sess__sv_%s_whiskerNoTouchTrials.mat',anim, date_str, sv_str_s2)); 
%          end 
         
         
       %% 
    %   start by finding MW Cells in one animal. day 1
         non_merged = find(~all_dat.anim_data{a}.is_merged_across_days);
         day_used_S1= non_merged(day_use);
         cts = {'usw', 'bsw', 'mw'};

            % restrict to desired area's id range
            %area_nrni = find(all_dat.anim_data{ai}.ids >= area_id_range{ai}(1) & all_dat.anim_data{ai}.ids <= area_id_range{ai}(2));

           if (stim_locked==0)  
            for c=1:length(cts) 
                ct = cts{c};
                cti = intersect(s1i,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1});
                idi = find(ismember(s1.dat.cellIds, all_dat.anim_data{a}.ids(cti)));
                cti_corrs= s1.dat.corrMat(idi,idi); %gives corr matrix for touch type
                cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                mean_corr_AT(ai,c,1)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type, s1
                
                cti_s2= (intersect(s2i,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}));%-(min(s2i-1));
                idi = find(ismember(s2.dat.cellIds, all_dat.anim_data{a}.ids(cti_s2)));
                cti_corrs_s2= (s2.dat.corrMat(idi,idi));
                cti_corrs_reflect_s2= triu(cti_corrs_s2)+triu(cti_corrs_s2,1)';
                cti_corrs_per_cell_s2= nanmean(cti_corrs_reflect_s2, 1);
                
                mean_corr_AT(ai,c,2)= mean(cti_corrs_per_cell_s2) ;
            end
            
            s1ci = find(ismember(s1.dat.cellIds, all_dat.anim_data{a}.ids(s1i)));
            mean_corr_AT(ai,4,1) = nanmean(reshape(s1.dat.corrMat(s1ci,s1ci),[],1)) ; % ALL cells
             
            s2ci = find(ismember(s2.dat.cellIds, all_dat.anim_data{a}.ids(s2i)));
            mean_corr_AT(ai,4,2) = nanmean(reshape(s2.dat.corrMat(s2ci,s2ci),[],1)) ; % ALL cells
           end  
           
           if (stim_locked)
               for aa=1:2 
                   if aa==1 
                       si=s1i; subvol= sv_str;
                   else si=s2i; subvol= sv_str_s2;
                   end 
          
                for c=1:length(cts) 
                ct = cts{c};
                cti = intersect(si,all_dat.types_by_idx{a}.([ct '_by_day']){day_used_S1}); %gives cells in area ID range that are mW, now need to break up further by touch type
              

                for tt=1:4 % loop through 4 touch types: w1p, w1r, w2p, w2r. trying to find s1, mw cells that specifically respond to each touch type 
                    if tt==1 touch_str= 'whiskerW1PExcTouchTrials'; type='ev_touch_w1p_by_day'; end 
                    if tt==2 touch_str= 'whiskerW2PExcTouchTrials'; type='ev_touch_w2p_by_day'; end 
                    if tt==3 touch_str= 'whiskerW1RExcTouchTrials'; type='ev_touch_w1r_by_day'; end 
                    if tt==4 touch_str= 'whiskerW2RExcTouchTrials'; type='ev_touch_w2r_by_day';end 
               
                    
                    ia= load (sprintf('%s_2022_%s_sess__sv_%s_%s.mat',anim, date_str, subvol, touch_str));
                    resp_cells=intersect(cti, (intersect(all_dat.types_by_idx{a}.([type]){day_used_S1}, si)));
                    idi= find(ismember(ia.dat.cellIds, all_dat.anim_data{a}.ids(resp_cells))); %so these are the IDS in the corr mat of MW cells that respond to W1P! 
                    num_cell (c, tt, ai)= length(idi); %sanity check if we need cell # cutoffs
                    cti_corrs= ia.dat.corrMat(idi,idi); %gives corr matrix for touch type --> these are correlations with OTHER cells of the exact same type
                    cti_corrs_reflect= triu(cti_corrs)+triu(cti_corrs,1)'; %mirrors corr mat about the diagonal to more easily compute
                    cti_corrs_per_cell= nanmean(cti_corrs_reflect, 1);
                    mean_corr_AT(c,tt, ai)= mean(cti_corrs_per_cell); % mean correlation for each animal for each cell type and for each touch type in area of interest
                    mean_corr_ancurr (c,aa)= mean(cti_corrs_per_cell);
                end 
                
                
                end 
               if aa==1 S1_means=cat(2,mean_corr_ancurr(:,:),mean_corr_ancurr(:,:)); mean_plot_rlly_S1(:,ai)= nanmean(S1_means,2); end 
               if aa==2 S2_means=cat(2,mean_corr_ancurr(:,:),mean_corr_ancurr(:,:)); mean_plot_rlly_S2(:,ai)= nanmean(S2_means,2); end 
               end 
            
            %means= cat(2,mean_corr_ancurr(:,:,1),mean_corr_ancurr(:,:,2));
            %mean_plot_rlly(:,ai)= nanmean(means,2); %columns are animals, rows are cell type: usw/bsw/mw
            %mean_plot_rlly(:,ai)= nanmean(mean,2);
           end 
           end 
    
%%
    if (stim_locked==0)
     x1= [1.5 3.5 5.5 7.5]; x2= [2.5 4.5 6.5 8.5];   
     color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color ; 0.75 0.75 0.75];
     names={'usw', 'bsw', 'mw','nt'};
    for cti=1:4 
        
         %title (ax_mean, 'mean response for cells of each type: averaged across all types')
         bar (ax, x1(cti),nanmean(nonzeros(mean_corr_AT(:,cti,1))),'FaceColor',color_arr(cti,:)); hold on; %s1
         plot (ax, x1(cti),(nonzeros(mean_corr_AT(:,cti,1))),'ko')
         bar (ax, x2(cti),nanmean(nonzeros(mean_corr_AT(:,cti,2))),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
         plot (ax, x2(cti),(nonzeros(mean_corr_AT(:,cti,2))),'ko')
         ax.TickDir='out';ax.Box='off';
         
       A =(mean_corr_AT(:,cti,1));
       B = (mean_corr_AT(:,cti,2));
          [h,pv]=ttest(A,B);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A), nanmean(B), nanstd(B), length(find(~isnan(A))), pv));
       if cti==1
              C= nonzeros(mean_corr_AT(:,cti+1,1)); %s1
              D= nonzeros(mean_corr_AT(:,cti+2,1));%s1
              F= nonzeros(mean_corr_AT(:,cti+1,2));%s2
              G= nonzeros(mean_corr_AT(:,cti+2,2));%s2
              Z= nonzeros(mean_corr_AT(:,cti+3,2));%s2
              Y= nonzeros(mean_corr_AT(:,cti+3,1));%s1
      
               
          [h,pv]=ttest(A,C);
          pval(cti)= pv;
          disp(sprintf('. %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(C), nanstd(C), length(find(~isnan(A))), pv));
      
          [h,pv]=ttest(A,D);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+2}, nanmean(D), nanstd(D), length(find(~isnan(A))), pv));
           
          [h,pv]=ttest(A,Y);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+2}, nanmean(Y), nanstd(Y), length(find(~isnan(A))), pv));
                
          [h,pv]=ttest(B,F);
          pval(cti)= pv;
          disp(sprintf(' %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(F), nanstd(F), length(find(~isnan(B))), pv));
      
          [h,pv]=ttest(B,G);
          pval(cti)= pv;
          disp(sprintf(' %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+2}, nanmean(G), nanstd(G), length(find(~isnan(B))), pv));
          
          [h,pv]=ttest(B,Z);
          pval(cti)= pv;
          disp(sprintf(' %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+2}, nanmean(Z), nanstd(Z), length(find(~isnan(B))), pv));
                      
          end
          if cti==2
          
          J= nonzeros(mean_corr_AT(:,cti+1,1));%s1
          JJ= nonzeros(mean_corr_AT(:,cti+1,2));%s2
          K= nonzeros(mean_corr_AT(:,cti+2,1));%s1
          KK= nonzeros(mean_corr_AT(:,cti+2,2));%s2
              
          [h,pv]=ttest(A,J);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(J), nanstd(J), length(find(~isnan(A))), pv));
          
          [h,pv]=ttest(A,K);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(K), nanstd(K), length(find(~isnan(A))), pv));
          
                
          [h,pv]=ttest(B,JJ);
          pval(cti)= pv;
          disp(sprintf('. %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(JJ), nanstd(JJ), length(find(~isnan(B))), pv));
          
          [h,pv]=ttest(B,KK);
          pval(cti)= pv;
          disp(sprintf('. %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(KK), nanstd(KK), length(find(~isnan(B))), pv));
                
          end 
          
          if cti==3
          
          V= nonzeros(mean_corr_AT(:,cti+1,1));%s1
          VV= nonzeros(mean_corr_AT(:,cti+1,2));%s2
          
          [h,pv]=ttest(A,V);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(V), nanstd(V), length(find(~isnan(A))), pv));
          
                
          [h,pv]=ttest(B,VV);
          pval(cti)= pv;
          disp(sprintf('. %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(VV), nanstd(VV), length(find(~isnan(B))), pv));
          
          end 
              
    end 

    end 
    if (stim_locked)
     x1= [1.5 3.5 5.5 ]; x2= [2.5 4.5 6.5 ];   
     color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
     names={'usw', 'bsw', 'mw'};
    for cti=1:3
        
         %title (ax_mean, 'mean response for cells of each type: averaged across all types')
         bar (ax, x1(cti),nanmean(nonzeros(mean_plot_rlly_S1(cti,:))),'FaceColor',color_arr(cti,:)); hold on; %s1
         plot (ax, x1(cti),(nonzeros(mean_plot_rlly_S1(cti,:))),'ko')
         bar (ax, x2(cti),nanmean(mean_plot_rlly_S2(cti,:)),'FaceColor',color_arr(cti,:),'FaceAlpha',0.4); hold on; %s1
         plot (ax, x2(cti),(nonzeros(mean_plot_rlly_S2(cti,:))),'ko')
         ax.TickDir='out';ax.Box='off';
         
       A = mean_plot_rlly_S1 (cti,:);
       B = mean_plot_rlly_S2 (cti,:);
          [h,pv]=ttest(A,B);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A), nanmean(B), nanstd(B), length(find(~isnan(A))), pv));
       if cti==1
              C= (mean_plot_rlly_S1 (cti+1,:)); %s1
              D= (mean_plot_rlly_S1 (cti+2,:));%s1
              F= (mean_plot_rlly_S2 (cti+1,:));%s2
              G= (mean_plot_rlly_S2 (cti+2,:));%s2
              %Z= nonzeros(mean_corr_AT(:,cti+3,2));%s2
              %Y= nonzeros(mean_corr_AT(:,cti+3,1));%s1
      
               
          [h,pv]=ttest(A,C);
          pval(cti)= pv;
          disp(sprintf('. %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(C), nanstd(C), length(find(~isnan(A))), pv));
      
          [h,pv]=ttest(A,D);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+2}, nanmean(D), nanstd(D), length(find(~isnan(A))), pv));
           
%           [h,pv]=ttest(A,Y);
%           pval(cti)= pv;
%           disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
%                     names{cti}, nanmean(A), nanstd(A),names{cti+2}, nanmean(Y), nanstd(Y), length(find(~isnan(A))), pv));
                
          [h,pv]=ttest(B,F);
          pval(cti)= pv;
          disp(sprintf(' %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(F), nanstd(F), length(find(~isnan(B))), pv));
      
          [h,pv]=ttest(B,G);
          pval(cti)= pv;
          disp(sprintf(' %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+2}, nanmean(G), nanstd(G), length(find(~isnan(B))), pv));
          
%           [h,pv]=ttest(B,Z);
%           pval(cti)= pv;
%           disp(sprintf(' %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
%                     names{cti}, nanmean(B), nanstd(B),names{cti+2}, nanmean(Z), nanstd(Z), length(find(~isnan(B))), pv));
%                       
          end
          if cti==2
          
          J= (mean_plot_rlly_S1(cti+1,:));%s1
          JJ= (mean_plot_rlly_S2(cti+1,:));%s2
          %K= nonzeros(mean_corr_AT(:,cti+2,1));%s1
          %KK= nonzeros(mean_corr_AT(:,cti+2,2));%s2
              
          [h,pv]=ttest(A,J);
          pval(cti)= pv;
          disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(J), nanstd(J), length(find(~isnan(A))), pv));
          
%           [h,pv]=ttest(A,K);
%           pval(cti)= pv;
%           disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
%                     names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(K), nanstd(K), length(find(~isnan(A))), pv));
%           
                
          [h,pv]=ttest(B,JJ);
          pval(cti)= pv;
          disp(sprintf('. %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
                    names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(JJ), nanstd(JJ), length(find(~isnan(B))), pv));
          
%           [h,pv]=ttest(B,KK);
%           pval(cti)= pv;
%           disp(sprintf('. %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
%                     names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(KK), nanstd(KK), length(find(~isnan(B))), pv));
%                 
          end 
          
%           if cti==3
%           
%           V= nonzeros(mean_corr_AT(:,cti+1,1));%s1
%           VV= nonzeros(mean_corr_AT(:,cti+1,2));%s2
%           
%           [h,pv]=ttest(A,V);
%           pval(cti)= pv;
%           disp(sprintf(' %s S1 mean/sd: %0.3f/%0.3f %s S1: %0.3f/%0.3f n=%d pval: %0.3f', ...
%                     names{cti}, nanmean(A), nanstd(A),names{cti+1}, nanmean(V), nanstd(V), length(find(~isnan(A))), pv));
%           
%                 
%           [h,pv]=ttest(B,VV);
%           pval(cti)= pv;
%           disp(sprintf('. %s S2 mean/sd: %0.3f/%0.3f %s S2: %0.3f/%0.3f n=%d pval: %0.3f', ...
%                     names{cti}, nanmean(B), nanstd(B),names{cti+1}, nanmean(VV), nanstd(VV), length(find(~isnan(B))), pv));
%           
%           end 
              
    end 
        
    end 
end
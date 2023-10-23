%  LR: 1/31/23: Figure that shows how much touch information is carried by
%  each set of projecting neurons of a given touch type: usw, bsw, mw. 

    all_dat = get_s1s2_all_dat;
    
    fh = figure('Position', [0 0 2000 500]); 
%     
    s1_to_s2_ax = subplot('Position', [.05 .05 .15 .8]); %set up plots 
    hold on;
    s2_to_s1_ax = subplot('Position', [.25 .05 .15 .8]);  
    hold on;

    vali = find(all_dat.retrograde_label_type ==  1); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    plot_proj_single_set(all_dat, vali, all_dat.s2_id_range, s1_to_s2_ax, 's1p neurons in s2');

    vali = find(all_dat.retrograde_label_type ==  2); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    plot_proj_single_set(all_dat, vali, all_dat.s1_id_range, s2_to_s1_ax,'s2p neurons in s1');
%     

  function plot_proj_single_set(all_dat, vali, area_id_range,ax,tstr);
    cts = {'usw', 'bsw', 'mw'}; %nt is no touch
    
    for v= 1:1:length(vali)
        ai=vali(v); 

        % restrict to desired area's id range
        nrni = find(all_dat.anim_data{ai}.ids >= area_id_range{ai}(1) & all_dat.anim_data{ai}.ids <= area_id_range{ai}(2));

        % map dates onto indexing
        d_num = nan*zeros(1,length(all_dat.anim_data{ai}.date_str));
        for d=1:length(all_dat.anim_data{ai}.date_str)
            if (length(all_dat.anim_data{ai}.date_str{d}) == 5)
                d_num(d) = find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}));
            end
        end
        di = min(find(all_dat.retro_dates{ai}));
        dat_dayi = find(d_num == di);
        
        %find which cells are responsive to which type of touch, make
        %matrix
        w1p_idx= intersect(all_dat.types_by_idx{ai}.ev_touch_w1p_by_day{dat_dayi}, nrni); %find indexes of each type of touch in s1 and s2
        w1r_idx= intersect(all_dat.types_by_idx{ai}.ev_touch_w1r_by_day{dat_dayi}, nrni);
        w2p_idx= intersect(all_dat.types_by_idx{ai}.ev_touch_w2p_by_day{dat_dayi}, nrni);
        w2r_idx= intersect(all_dat.types_by_idx{ai}.ev_touch_w2r_by_day{dat_dayi}, nrni);    
       
        for c=1:length(cts) ; %go trhough and perform computations for each group of cells: usw, bsw, mw.
            ct = cts{c};
            cti = intersect(nrni,all_dat.types_by_idx{ai}.([ct '_by_day']){dat_dayi}); %gives cells that are in the id range and of a given type for the given dates
            N_red_ct(c) = length(intersect(cti,all_dat.types_by_idx{ai}.is_red));
            N_ct(c) = length(cti);
            red_cells_ids= intersect(cti,all_dat.types_by_idx{ai}.is_red) %red cells in a given group on a given day 

            %red_total= length(all_dat.types_by_idx{ai}.is_red);
%             
            resp_mean= nan*zeros(length(red_cells_ids),4);
            
            touch_types =nan*zeros(length(red_cells_ids), 4);
            for i= 1: length(red_cells_ids) %loop thorugh the possible all cell ids and make a matrix of cell id x type (w1p, w1r, w2p, w2r) per animal in S1 
               touch_types(i,1)= any(w1p_idx(:)== red_cells_ids(i));
               touch_types(i,2)= any(w1r_idx(:)== red_cells_ids(i));
               touch_types(i,3)= any(w2p_idx(:)== red_cells_ids(i));
               touch_types(i,4)= any(w2r_idx(:)== red_cells_ids(i));   
            end
            
            for tt= 1:4 % loop through 4 touch types: w1p, w1r, w2p, w2r
               if tt==1 response=all_dat.types_by_idx{ai}.ev_touch_w1p_by_day, mean_use=all_dat.anim_mats{ai}.meanRespW1p; prob_use=all_dat.anim_mats{ai}.probRespW1p; end 
               if tt==3 response=all_dat.types_by_idx{ai}.ev_touch_w2p_by_day; mean_use= all_dat.anim_mats{ai}.meanRespW2p; prob_use=all_dat.anim_mats{ai}.probRespW2p;end 
               if tt==2 response=all_dat.types_by_idx{ai}.ev_touch_w1r_by_day; mean_use= all_dat.anim_mats{ai}.meanRespW1r;prob_use=all_dat.anim_mats{ai}.probRespW1r;end 
               if tt==4 response=all_dat.types_by_idx{ai}.ev_touch_w2r_by_day; mean_use= all_dat.anim_mats{ai}.meanRespW2r;prob_use=all_dat.anim_mats{ai}.probRespW2r;end 
               if N_red_ct(c)==0 mean_use=zeros(1,5000); red_cells_ids=1; touch_types= zeros(1,4); resp_mean=0; end
            

              resp_mean(:, tt) = mean_use(dat_dayi, red_cells_ids)'; %this is the means of red cell ids for that touch type
          
           end 
           
           resp_mean(isnan(resp_mean))=1; % need someway to non let nans absolutely derail everything--doesnt effect overall computation 
           
           for jj= 1: length(red_cells_ids)
               mean_per_cell(jj)= mean(nonzeros(resp_mean(jj,:).*touch_types(jj,:))); %find mean response of each cell to each touch type  
           end 
          
           sum_report (ai, c)= sum(mean_per_cell) %gives sum of the mean of each projecting cell for each category usw/bsw/mw for each animal
       
        end
        
        N_touch = sum(N_ct);
        N_all = length(nrni);
        N_red = length(all_dat.types_by_idx{ai}.is_red);
        
        ratio_report(ai,:)= N_red_ct./N_red; %this is fraction of projecting cells of that type out of all projecting cells
        %ratio_report(ai,4)= (N_red-sum(N_red_ct(1:3)))/N_red; %what
        %fraction red cells is non touch, so the sum of ratios for each
        %animal should be 1 -->at this point, not sure what this adds
        norm_ratio (ai, :)= ratio_report(ai,:)/sum(ratio_report(ai,:)); %not using these ratios right now, might use eventually
       
    end 
    
sum_report(isnan(sum_report))=0; %make nans 0 since they contribute 0 to overall response
overall_sum= sum(sum_report,2); %find sum of ALL touch activity in the area
sum_normed= sum_report./overall_sum; %this should give what % of touch activity is carried by which cells

x1= [1.5 3.5 5.5]; %plot 
color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color];
% 
   for cts=1:3   
       bar (ax, x1(cts),mean((sum_normed(:,cts)),'omitnan'),'FaceColor',color_arr(cts,:)); hold on; %s1
       plot (ax, x1(cts),((sum_normed(:,cts))),'ko')
       title (ax, tstr); 
       set(ax,'TickDir','out','FontSize',15);
   end 
 end
 
          
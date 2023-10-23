% LR updated 9/21/23
% code to find identity shifts in baseline two days and also pre/post
% lesion

all_dat = get_s1s2_all_dat; %setup
%% USE FOR S2 lesion/sham
 
%   s2_lesion_ani = find(~strcmp('',all_dat.s2_lesion_date));
%   s2_lesion_ani = s2_lesion_ani(2:end); %animal id 12 is messed up for the baseline (dont have enough data to compare days)
%   frac_s2_lesion = plot_identity_shift_lesion(all_dat, s2_lesion_ani, all_dat.lesion_dates, all_dat.s1_lesion_id_range, 's2 lesion');
%   
%   sham_lesion_ani = find(all_dat.is_1daysham); % restrict to only 1day shams (This is more appropriate - what we do for everything else!)
%   frac_sham= plot_identity_shift_lesion(all_dat, sham_lesion_ani, all_dat.sham_dates, all_dat.s1_lesion_id_range, 'Sham lesion');      
% 
%   s1_baseline_ani= s2_lesion_ani; 
%   frac_s1_baseline= plot_identity_shift_lesion(all_dat, s1_baseline_ani, '', all_dat.s1_lesion_id_range, 'Baseline s1')
%     
%  % once function runs for these and we get plots, then we need to calculate
%  % all the p values 
%     for i= 1:4
%       for j=1:4
%          [~,pv(i,j)]= ttest(frac_s1_baseline(i,j,:), frac_s2_lesion(i,j,:)) %paired
%          [~,psham(i,j)]= ttest2(frac_s1_baseline(i,j,:), frac_sham(i,j,:)) %unpaired
%       end 
%     end 

 %% USE section for S1 lesion

 s1_lesion_ani = find(~strcmp('',all_dat.s1_lesion_date));     
 s1_lesion_ani= s1_lesion_ani([3:end]); %17518 has too few touch cells for thsi
 frac_s1_lesion= plot_identity_shift_lesion(all_dat, s1_lesion_ani, all_dat.lesion_dates, all_dat.s2_lesion_id_range, 's1 lesion');

 s2_baseline_ani= s1_lesion_ani; 
 frac_s2_baseline= plot_identity_shift_lesion(all_dat, s2_baseline_ani, '', all_dat.s2_lesion_id_range, 'Baseline s2')
  
  for i= 1:4
      for j=1:4
         [~,pv(i,j)]= ttest(frac_s2_baseline(i,j,:), frac_s1_lesion(i,j,:)) %paired
      end 
  end 
%% Function to calculate and plot the identity shifts, outputs frac of each pre day identitiy in each post day identity group 
  
function frac_pre_post= plot_identity_shift_lesion(all_dat, ani_used, dates_used,id_range,les_type)
    
    for ai=1:length(ani_used)
        a = ani_used(ai);
        ndays = length(all_dat.anim_data{a}.date_str);

        % restrict by ID 
        ni = find(all_dat.anim_data{a}.ids >= id_range{a}(1) & all_dat.anim_data{a}.ids <= id_range{a}(2));
        if (length(ni) == 0) ; disp([' **** PROBLEM **** ' all_dat.anims{a} ' has ZERO neurons ']) ; end
        
        % which date should we use?
        if strcmp(les_type, 'Baseline s2') || strcmp(les_type, 'Baseline s1')
            usw_id_pre = intersect(ni,all_dat.types_by_idx{a}.usw_by_day{1}); %usw cell ids prior to lesion
            bsw_id_pre = intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{1});
            mw_id_pre= intersect(ni,all_dat.types_by_idx{a}.mw_by_day{1});
            nt_id_pre= setdiff(setdiff((setdiff(ni, usw_id_pre)), bsw_id_pre),mw_id_pre);% no touch cells ids prior to leison
            usw_id_post = intersect(ni,all_dat.types_by_idx{a}.usw_by_day{2});
            bsw_id_post = intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{2});
            mw_id_post = intersect(ni,all_dat.types_by_idx{a}.mw_by_day{2}); 
        else
            pre_di = find(dates_used{a} == -1);
            post_di = find(dates_used{a} == 1);

            % map onto indexing
            d_num = nan*zeros(1,length(all_dat.anim_data{a}.date_str));
            for d=1:length(all_dat.anim_data{a}.date_str)
                if (length(all_dat.anim_data{a}.date_str{d}) == 5)
                    d_num(d) = find(strcmp(all_dat.valid_dates{a}, all_dat.anim_data{a}.date_str{d}));
                end
            end
        
        
        for d=1:length(pre_di)
            di = find(d_num == pre_di(d));
            usw_id_pre = intersect(ni,all_dat.types_by_idx{a}.usw_by_day{di}); %usw cell ids prior to lesion
            bsw_id_pre = intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{di});
            mw_id_pre= intersect(ni,all_dat.types_by_idx{a}.mw_by_day{di});
            nt_id_pre= setdiff(setdiff((setdiff(ni, usw_id_pre)), bsw_id_pre),mw_id_pre);% no touch cells ids prior to leison
        end
         
        
        for d=1:length(post_di)
            di = find(d_num == post_di(d));
            usw_id_post = intersect(ni,all_dat.types_by_idx{a}.usw_by_day{di});
            bsw_id_post = intersect(ni,all_dat.types_by_idx{a}.bsw_by_day{di});
            mw_id_post = intersect(ni,all_dat.types_by_idx{a}.mw_by_day{di}); 
        end
        
        end
        
        
        
        sum_cells= length(ni);
        touch_pre= length(usw_id_pre) +length(bsw_id_pre) +length(mw_id_pre);
        touch_post= length(usw_id_post) +length(bsw_id_post) +length(mw_id_post);
        type_nums= [length(usw_id_pre) length(bsw_id_pre) length(mw_id_pre) length(nt_id_pre)];
        frac_per_an(ai,:)= type_nums./sum_cells;
        
        cts = {'usw', 'bsw', 'mw','nt'}; %nt is no touch
        for c=1:length(cts) ; %go trhough and perform computations for each group of cells: usw, bsw, mw.
            ct = cts{c};
            
            if c==1 pre_id= usw_id_pre; end 
            if c==2 pre_id= bsw_id_pre; end
            if c==3 pre_id= mw_id_pre;  end
            if c==4 pre_id= nt_id_pre;  end 
            
            
            id_prepost(c,1,ai)= length(intersect(pre_id, usw_id_post)); %this is the number of cells of type c that switch to usw after lesion
            id_prepost(c,2,ai)= length(intersect(pre_id, bsw_id_post)); %num cells that start at type c that become type bsw after lesion
            id_prepost(c,3,ai)= length(intersect(pre_id, mw_id_post)); %num cells type c that become mw after lesion
            id_prepost(c,4,ai)= length(setdiff((setdiff((setdiff(pre_id,(intersect(pre_id, usw_id_post)))),(intersect(pre_id, bsw_id_post)))),(intersect(pre_id, mw_id_post))));%num cells type c that beocme non touch after lesion
            
          
        end 
        
        frac_pre_post(:,:,ai)= id_prepost(:,:,ai)./sum(id_prepost(:,:,ai),2); %individual animals, turnover fraction
%         frac_touch_pre(ai,:)= [length(usw_id_pre)/touch_pre length(bsw_id_pre)/touch_pre length(mw_id_pre)/touch_pre];
%         frac_touch_post(ai,:)= [length(usw_id_post)/touch_post length(bsw_id_post)/touch_post length(mw_id_post)/touch_post];
        frac_touch_pre(ai,:)= [length(usw_id_pre)/sum_cells length(bsw_id_pre)/sum_cells length(mw_id_pre)/sum_cells];
        frac_touch_post(ai,:)= [length(usw_id_post)/sum_cells length(bsw_id_post)/sum_cells length(mw_id_post)/sum_cells];
 
        for cc=1:4
            norm_test(cc,:,ai)=  frac_pre_post(cc,:,ai).*frac_per_an(ai,cc);
        end 
    end
    
    mean_across_animals_frac= mean(frac_pre_post,3);
    mean_norm_test= mean(norm_test,3); 
    mean_across_an_pre= mean(frac_touch_pre,1); %mean proportion of touch pop that usw/bsw/mw make up across animals
    mean_across_an_post= mean(frac_touch_post,1);
    
%     
     x1= [1.5 3.5 5.5]; %plot 
        color_arr= [all_dat.usw_color; all_dat.bsw_color; all_dat.mw_color; [0 0 0]];
       
        fx = figure('Position', [0 0 500 500]); 
        %ax2=bar(mean_across_animals_frac,'stacked')
        ax2= bar(mean_across_animals_frac(:,:), 'stacked')
        ax2(1).FaceColor= color_arr(1,:);
        ax2(2).FaceColor= color_arr(2,:);
        ax2(3).FaceColor= color_arr(3,:);
        ax2(4).FaceColor= color_arr(4,:);
        
        xlabel ('type pre lesion')
        xticklabels({'usw', 'bsw', 'mw', 'nt'})
        title (les_type)
        
end 
        
    
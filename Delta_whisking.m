%LR 5/24/23, whisking test

all_dat = get_s1s2_all_dat;
s1_lesion_ani = find(~strcmp('',all_dat.s1_lesion_date));
s2_lesion_ani = find(~strcmp('',all_dat.s2_lesion_date));




%%
[delt_whisk_s1les, pv_s1les, anid_s1les, f_pre_s1les, f_post_s1les]= delta_whisk(1, all_dat, s1_lesion_ani, 's2_manual_restrict', all_dat.lesion_dates, 's1_lesion');

[delt_whisk_s2les, pv_s2les, anid_s2les, f_pre_s2les, f_post_s2les]= delta_whisk(2, all_dat, s2_lesion_ani, 's1_manual_restrict', all_dat.lesion_dates,'s2_lesion'); 
 %% for the moment 
 %ani= s1_lesion_ani; 
 %restrict_tag= 's2_manual_restrict';
 %dates_used= all_dat.lesion_dates; 
 %tstr= 's1_lesion';
 %i=1;
%% 
function [delt_whisk, pv, delt_ani, frac_whisk_pre,frac_whisk_post]= delta_whisk(i, all_dat, ani, restrict_tag, dates_used,tstr)
    for ai= 1: length(ani)
        a = ani(ai);
        ndays = length(all_dat.anim_data{a}.date_str);

            % restrict by ID 
            ni = get_s1s2_neuron_subset_idx(all_dat.anim_data{a}.ids, restrict_tag, all_dat, a);


            % which date should we use?
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
                whisk_cells= intersect(ni, all_dat.types_by_idx{a}.ev_whisking_by_day{di});
                pre_whisk(d)= mean(all_dat.anim_mats{a}.meanRespWhisk(di,whisk_cells));
            end 

            pre_respw(ai)= mean(pre_whisk); 
            frac_whisk_pre(ai)= length(whisk_cells)/length(ni);


            for d=1:length(post_di)
                di = find(d_num == post_di(d));
                whisk_cells= intersect(ni, all_dat.types_by_idx{a}.ev_whisking_by_day{di});
                post_whisk(d)= mean(all_dat.anim_mats{a}.meanRespWhisk(di,whisk_cells));
            end 

            post_respw(ai)= mean(post_whisk); 
            frac_whisk_post(ai)= length(whisk_cells)/length(ni);
            
           

    end 
    delt_whisk= mean(post_respw- pre_respw);
    delt_ani= (post_respw- pre_respw);
    %normd_delt_ani= (post_respw- pre_respw)./pre_respw;
    %normd_mean= mean(normd_delt_ani);
    [h, pv]= ttest(pre_respw-post_respw); 
   
    %plot(ax,[0.5 0.6], [1 1]*delt_whisk, 'k-', 'LineWidth', 5);
    %plot(ax,0.55*ones(1,length(ani)), delt_ani, 'o', 'Color', [1 1 1]*0.25, 'MarkerFaceColor', [1 1 1]*0.75);
   
    %plot (ax, [1 2], [mean(pre_respw), mean(post_respw)], 'k-', 'LineWidth',5); hold on;
    %for ii=1:length(ani)
        %plot (ax, [1 2], [pre_respw(ii), post_respw(ii)], 'm--');
   % end
    %
    %ylim(ax, [-0.5 1]);
    %xlim(ax, [0 3]);
    %t%itle(ax,tstr);
    x= [1 2 3 4];
    figure (i);
    ax1= subplot (1,2,1); hold on;
    ax2= subplot (1,2,2); hold on;
    bar (ax1, x(1:2), [mean(frac_whisk_pre),mean(frac_whisk_post)]); hold on;
    for ii=1:length(ani)
         plot (ax1, [1 2], [frac_whisk_pre(ii), frac_whisk_post(ii)], 'ko');
    end
    this_tstr= ['fraction whisk', tstr];
    title (ax1,this_tstr);
    %ylim(ax1, [0 0.5])
    
    bar (ax2, x(1:2), [mean(pre_respw), mean(post_respw)]); 
    for ii=1:length(ani)
            plot (ax2, [1,2], [pre_respw(ii), post_respw(ii)], 'ko'); 
    end 
    this_tstr= ['mean whisk', tstr];
    title (ax2,this_tstr)
    %ylim(ax2, [0 0.5])

end


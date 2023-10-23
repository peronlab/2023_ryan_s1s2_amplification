function figure_s1s2_corr
    fh = figure('Position',[0 0 900 300]);
    s1_corr_ax = subplot('Position', [.1 .1 .25 .75]);
    s2_corr_ax = subplot('Position', [.4 .1 .25 .75]);

    all_dat = get_s1s2_all_dat;;

    anis = find(all_dat.has_deep_data_s1s2 ==  1); % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    plot_corr_single_set(all_dat, anis, all_dat.s1_id_range, all_dat.static_dates, 1, s1_corr_ax, 's1 dual volume');
    plot_corr_single_set(all_dat, anis, all_dat.s2_id_range, all_dat.static_dates, 1, s2_corr_ax, 's2 dual volume');

    anis = find(~strcmp(all_dat.s1_lesion_date, ''));
    plot_corr_single_set(all_dat, anis, all_dat.s2_lesion_id_range, all_dat.lesion_dates, -1, s1_corr_ax, 's1 lesion pre');
    plot_corr_single_set(all_dat, anis, all_dat.s2_lesion_id_range, all_dat.lesion_dates, 1, s2_corr_ax, 's1 lesion post');

    anis = find(~strcmp(all_dat.s2_lesion_date, ''));
    plot_corr_single_set(all_dat, anis, all_dat.s1_lesion_id_range, all_dat.lesion_dates, -1, s1_corr_ax, 's2 lesion pre');
    plot_corr_single_set(all_dat, anis, all_dat.s1_lesion_id_range, all_dat.lesion_dates, 1, s2_corr_ax, 's2 lesion post');



%    axis(s1_to_s2_ax, [aa_12(1:3) max(aa_12(4),aa_21(4))]);
%    axis(s1_corr_ax, [aa_21(1:3) max(aa_12(4),aa_21(4))]);


function plot_corr_single_set(all_dat, anis, area_id_range, dates_passed, date_value_used, ax, tstr);
    cts = {'usw', 'bsw', 'mw'};
    
    disp(['--------------- ' tstr ' ----------------']);

    mu_corrs_mat = nan*zeros(3,length(anis), 2); % cell type X animal X alltime/notouch
    for a=1:length(anis)
        ai=anis(a); 
        nrni = find(all_dat.anim_data{ai}.ids >= area_id_range{ai}(1) & all_dat.anim_data{ai}.ids <= area_id_range{ai}(2));

        % map dates onto indexing ; also load corr data (!)
        d_num = nan*zeros(1,length(all_dat.anim_data{ai}.date_str));
        for d=1:length(all_dat.anim_data{ai}.date_str)
            if (length(all_dat.anim_data{ai}.date_str{d}) == 5)
                if(~isempty(find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}))));
                    d_num(d) = find(strcmp(all_dat.valid_dates{ai}, all_dat.anim_data{ai}.date_str{d}));
                end
            end
        end

        % only use specified date
        dates_used = dates_passed{ai};
        datei = find(dates_used == date_value_used);
        dat_dir = [all_dat.root_dir filesep all_dat.anims{ai} filesep 'session_neuropilone/pairwiseCorrelations'];
        submat = nan*zeros(3,2,length(datei));
        for d=1:length(datei)
            di = find(d_num == datei(d));

            % get data
            all_time_fs = dir([dat_dir filesep '/*' all_dat.anim_data{ai}.date_str{di} '*_allTime.mat']);
            if (isempty(all_time_fs)) ; disp(['NO MATCH: ' dat_dir filesep '/*' all_dat.anim_data{ai}.date_str{di} '*_allTime.mat']); end
            tmp_corr_mat = nan*zeros(length(all_dat.anim_data{ai}.ids));
            for f=1:length(all_time_fs)
            %    disp(['Loading: ' dat_dir filesep all_time_fs(f).name]);
                load([dat_dir filesep all_time_fs(f).name]);
                ii = find(ismember(all_dat.anim_data{ai}.ids, dat.cellIds));
                valci = find(ismember(dat.cellIds, all_dat.anim_data{ai}.ids(nrni)));
                if (length(ii) > 0 & length(valci) == length(dat.cellIds))
                    tmp_corr_mat(ii,ii) = dat.corrMat;
                end

                if(length(find(isnan(dat.corrMat(:)))) == length(dat.corrMat(:))) ; disp(['nan BAD: ' dat_dir filesep all_time_fs(f).name]); end
            end

            % select cells
            uswi = all_dat.types_by_idx{ai}.usw_by_day{di};
            bswi = all_dat.types_by_idx{ai}.bsw_by_day{di};
            mwi = all_dat.types_by_idx{ai}.mw_by_day{di};

            submat(1,1,d) = nanmean(nanmean(tmp_corr_mat(uswi,uswi)));
            submat(2,1,d) = nanmean(nanmean(tmp_corr_mat(bswi,bswi)));
            submat(3,1,d) = nanmean(nanmean(tmp_corr_mat(mwi,mwi)));

            notouch_time_fs = dir([dat_dir filesep '/*' all_dat.anim_data{ai}.date_str{di} '*whiskerNoTouchTrials.mat']);
            if (isempty(notouch_time_fs)) ; disp(['NO MATCH: ' dat_dir filesep '/*' all_dat.anim_data{ai}.date_str{di} '*whiskerNoTouchTrials.mat']); end
            tmp_corr_mat = nan*tmp_corr_mat;
            for f=1:length(notouch_time_fs)
            %    disp(['Loading: ' dat_dir filesep notouch_time_fs(f).name]);
                load([dat_dir filesep notouch_time_fs(f).name]);
                ii = find(ismember(all_dat.anim_data{ai}.ids, dat.cellIds));
                valci = find(ismember(dat.cellIds, all_dat.anim_data{ai}.ids(nrni)));
                if (length(ii) > 0 & length(valci) == length(dat.cellIds))
                    tmp_corr_mat(ii,ii) = dat.corrMat;
                end

                if(length(find(isnan(dat.corrMat(:)))) == length(dat.corrMat(:))) ; disp(['nan BAD: ' dat_dir filesep all_time_fs(f).name]); end
            end

            submat(1,2,d) = nanmean(nanmean(tmp_corr_mat(uswi,uswi)));
            submat(2,2,d) = nanmean(nanmean(tmp_corr_mat(bswi,bswi)));
            submat(3,2,d) = nanmean(nanmean(tmp_corr_mat(mwi,mwi)));
        end

        for c=1:3
            mu_corrs_mat(c,a,1) = nanmean(squeeze(submat(c,1,:)));
            mu_corrs_mat(c,a,2) = nanmean(squeeze(submat(c,2,:)));
        end

    end

mu_corrs_mat(:,:,1)'
nanmean(mu_corrs_mat(:,:,1)')

mu_corrs_mat(:,:,2)'
nanmean(mu_corrs_mat(:,:,2)')
        % restrict to desired area's id range

        % day to use -- not ideal
%        di = min(find(all_dat.retro_dates{ai}));
%        dat_dayi = find(d_num == di);
%        disp( ' ' ); 

    ms = 7;
    
    hold(ax, 'on');
    title(ax, tstr);
 if (0)   
    plot(ax,0*ratio_mat(1,:)+1, ratio_mat(1,:), 'o', 'Color', all_dat.usw_color, 'MarkerSize', ms);
    plot(ax,1, nanmean(ratio_mat(1,:)), 'o', 'Color', all_dat.usw_color, 'MarkerSize', ms+4, 'LineWidth',2);
    
    plot(ax,0*ratio_mat(2,:)+2, ratio_mat(2,:), 'o', 'Color', all_dat.bsw_color, 'MarkerSize', ms);
    plot(ax,2, nanmean(ratio_mat(2,:)), 'o', 'Color', all_dat.bsw_color, 'MarkerSize', ms+4, 'LineWidth',2);

    plot(ax,0*ratio_mat(3,:)+3, ratio_mat(3,:), 'o', 'Color', all_dat.mw_color, 'MarkerSize', ms);
    plot(ax,3, nanmean(ratio_mat(3,:)), 'o', 'Color', all_dat.mw_color, 'MarkerSize', ms+4, 'LineWidth',2);
    end
    set(ax,'TickDir','out','FontSize',15);

    aa=axis(ax);
    axis(ax, [0 4 0 aa(4)]);



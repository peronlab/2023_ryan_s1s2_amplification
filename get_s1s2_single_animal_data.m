%
% Gathers data for s1s2 project from session object -- only include stuff that has to access session objects here
%
function dat = get_s1s2_single_animal_data(force_redo)
    if (nargin == 0 ) ; force_redo = 0 ; end

    % preliminaries
    mdat = get_s1s2_metadata;
    fname = [pwd filesep 's1s2_summary_data.mat'];
    cwd_ani_str = fileparts_crossplatform(fileparts_crossplatform(pwd),2);
    anii = find(ismember(mdat.anims, cwd_ani_str));

    if (isempty(anii))
        disp('Current directory does not belong to an animal in metadata; add it to get_s1s2_metadata');
        return;
    end

    % main call
    if (exist(fname) & ~force_redo)
        load(fname);
    else
        if (force_redo)
            disp([fname ' getting regenerated - force_redo flag 1']);
        else
            disp([fname ' not found; regenerating from scratch']);
        end
        fl = dir('an*sess.mat');

        % pull ids from object #1, 2
        dat.ids = [];
        for f=1:length(fl) % we could go higher here if we find this is not enough ...
            % make sure we want this day
            uidx = find(fl(f).name == '_');
            date_mmdd_str = fl(f).name(uidx(2)+1:uidx(4)-1);

            if (ismember(date_mmdd_str, mdat.valid_dates{anii}))
                load(fl(f).name);
                if (isfield(s.caTSAArray,'caTSA'))
                    for c=1:length(s.caTSAArray.caTSA)
                        dat.ids = union(dat.ids, s.caTSAArray.caTSA{c}.ids);
                    end
                else
                    dat.ids = union(dat.ids, s.caTSA.ids);
                end
            else 
                disp(['For animal ' mdat.anims{anii} ', date ' date_mmdd_str ' is not part of valid_dates; skipping']);
            end
        end

        % preassign xyz
        dat.x_um = nan*dat.ids;
        dat.y_um = nan*dat.ids;
        dat.z_um = nan*dat.ids;
;

        % compoounds? these are always for early days so put them first. we take these if metadata says so
        f_base = 0;
        if (mdat.load_merged)
            fl = dir('../session_neuropilone_final_merged/an*sess.mat');
            if (~isempty(fl))
                opwd = pwd;
                cd (['../session_neuropilone_final_merged/']);
                for f=1:length(fl)
                    dat = process_single_sess(fl(f).name, f, dat, mdat, anii);
                end
                cd (opwd);
                f_base =f;
            end
        end

        % regular objects
        fl = dir('an*sess.mat');
        fi = 1;
        for f=1:length(fl)
            % make sure we want this day
            uidx = find(fl(f).name == '_');
            date_mmdd_str = fl(f).name(uidx(2)+1:uidx(4)-1);

            if (ismember(date_mmdd_str, mdat.valid_dates{anii})) % only take valid dates
                dat = process_single_sess(fl(f).name, fi+f_base, dat, mdat, anii);
                fi = fi+1;
            end
        end

        save(fname, 'dat');
    end

function dat = process_single_sess(fname, f, dat, mdat, anii)
    disp([ 'Processing  ' fname]);
    load(fname);
    if (isfield(s.caTSAArray,'caTSA'))
        this_ids = [];
        for c=1:length(s.caTSAArray.caTSA)
            this_ids = [this_ids s.caTSAArray.caTSA{c}.ids];
        end
    else
        try
            this_ids = s.caTSA.ids;
        catch me
            disp('SKIPPING!!!');
        end
    end
    olap_frac = (length(intersect(this_ids, dat.ids))/length(union(this_ids,dat.ids)));
    disp(['  ID overlap: ' num2str(olap_frac) ' (1 = all IDs present; <1, only subset)']);

    % evoked scores
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseTouchW1ExcludeW2');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseTouchW1ExcludeW2');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseTouchW2ExcludeW1');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseTouchW2ExcludeW1');

    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseTouchW1pExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseTouchW1pExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseTouchW1rExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseTouchW1rExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseTouchW2pExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseTouchW2pExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseTouchW2rExclusive');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseTouchW2rExclusive');
    
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_probResponseWhiskExcludeOptogenTouch');
    dat = add_single_evoked_feature(s, dat, f, 'evokedDffScores_meanResponseWhiskExcludeOptogenTouch');

    % touch model (add!)
    dat = get_two_whisker_model_fields (s, dat, f);

    % whisking GLM stuff
    dat = add_single_evoked_feature(s, dat, f, 'dffBased_k16t14_WhiskerSetpointExcludeTouchGLMScore');
    dat = add_single_evoked_feature(s, dat, f, 'dffBased_k16t14_MeanWhiskerAmplitudeExcludeTouchGLMScore');
    dat = add_single_evoked_feature(s, dat, f, 'dffBased_k16t14_WhiskerSetpointGLMScore');
    dat = add_single_evoked_feature(s, dat, f, 'dffBased_k16t14_MeanWhiskerAmplitudeGLMScore');

    % decoding (add!)

    dat.day(f).source_filename = fname;

    if (isfield(s.caTSAArray,'caTSA'))
        for c=1:length(s.caTSAArray.caTSA)

            % get x/y/z
            for fi=1:length(s.caTSAArray.caTSA{c}.roiArray)
                if (fi == 1)
                    % "base" z - topmost of this subvolume has a z of ???
                    %  1) are we s1/s2?
                    s1_or_s2 = 's2';
                    ri = find(dat.ids == s.caTSAArray.caTSA{c}.roiArray{fi}.roiIds(1));
                    if (dat.ids(ri) >= mdat.s1_id_range{anii}(1) & dat.ids(ri) <= mdat.s1_id_range{anii}(2) ) % s1!
                        s1_or_s2 = 's1';
                    end

                    % are we superficial/deep?
                    subvoli = str2num(s.caTSAArray.caTSA{c}.roiArray{fi}.metaData.subvolId(end));

                    % now we can infer which base z to use
                    if (strcmp(s1_or_s2,'s1') & subvoli == 1) ; 
                        base_z = mdat.s1_start_depth_super_rel_L1_micron(anii);
                    elseif (strcmp(s1_or_s2,'s1') & subvoli == 2) ; 
                        base_z = mdat.s1_start_depth_deep_rel_L1_micron(anii);
                    elseif (strcmp(s1_or_s2,'s2') & subvoli == 1) ; 
                        base_z = mdat.s2_start_depth_super_rel_L1_micron(anii);
                    elseif (strcmp(s1_or_s2,'s2') & subvoli == 2) 
                        base_z = mdat.s2_start_depth_deep_rel_L1_micron(anii);
                    else 
                        base_z = nan;
                    end
                    if (base_z == 0) ; base_z = nan; end % value of 0 is used for mice where this was not measured, so convert to nan
                end

                % COMs
                coms = s.caTSAArray.caTSA{c}.roiArray{fi}.getCOMs;
                %% roi loop w assignments
                for r=1:size(coms,1)
                    ri = find(dat.ids == s.caTSAArray.caTSA{c}.roiArray{fi}.roiIds(r));
                    if (length(ri) == 1 && isnan(dat.x_um(ri)))
                        dat.x_um(ri) = s.caTSAArray.caTSA{c}.roiArray{fi}.settings.umPerPixel(1)*coms(r,1);
                        dat.y_um(ri) = s.caTSAArray.caTSA{c}.roiArray{fi}.settings.umPerPixel(2)*coms(r,2);

                        % z is a bit different ...
                        y_frac = coms(r,2)/512;
                        dat.z_um(ri) = base_z + y_frac*mdat.dz_used(anii) + (fi-1)*mdat.dz_used(anii);
                    end
                end
            end 
        end   
    end
    

    % redness (for identifying retrograde cells)
    use_demix = 0; % seems raw works better
    use_lum_norm = 1; % luminance normalization before calculating redness?
    dat.day(f).red_score = nan*dat.ids;

    % loop over planes
    if (mdat.retrograde_label_type(anii) > 0)
        disp('Animal has retrograde label ; determining redness score.')
        if (isfield(s.caTSAArray,'caTSA'))
            for c=1:length(s.caTSAArray.caTSA)


                % demixed img
                dim = {[],[],[]};
                for ci=1:s.caTSAArray.caTSA{c}.length
                    fi = s.caTSAArray.caTSA{c}.roiFOVidx(ci);

                    % get demixed image
                    try
                        if (isempty(dim{fi}))
                            if (use_demix)
                                dim{fi} = s.caTSAArray.caTSA{c}.roiArray{fi}.demixRedChannel(s.caTSAArray.caTSA{c}.roiArray{fi}.masterImage);
                            elseif (use_lum_norm)
                                baseim = s.caTSAArray.caTSA{c}.roiArray{fi}.masterImage(:,:,2);
                                nim = normalize_via_gaussconvdiv(baseim, .025);
                                nim = nim/(quantile(nim(:),.99));
                                dim{fi} = nim;
                            else
                                dim{fi} = s.caTSAArray.caTSA{c}.roiArray{fi}.masterImage(:,:,2);
                            end
                        end

                        % score
                        roi = s.caTSAArray.caTSA{c}.roiArray{fi}.getRoiById(s.caTSAArray.caTSA{c}.ids(ci));
                        cci = find(dat.ids == s.caTSAArray.caTSA{c}.ids(ci));
                        if (length(cci) == 1 & ~isempty(roi))
                            dat.day(f).red_score(cci) = nanmean(dim{fi}(roi.indices));
                        end
                    catch me
                        disp('REDNESS FAILED!');
                    end
                end
            end
        else
            for ci=1:s.caTSA.length
                dim = {[],[],[]};
                fi = s.caTSA.roiFOVidx(ci);

                % get demixed image
                try
                    if (isempty(dim{fi}))
                        if (use_demix)
                            dim{fi} = s.caTSA.roiArray{fi}.demixRedChannel(s.caTSA.roiArray{fi}.masterImage);
                        elseif (use_lum_norm)
                            baseim = s.caTSA.roiArray{fi}.masterImage(:,:,2);
                            nim = normalize_via_gaussconvdiv(baseim, .025);
                            nim = nim/(quantile(nim(:),.99));
                            dim{fi} = nim; 
                        else % just normalize
                            dim{fi} = s.caTSA.roiArray{fi}.masterImage(:,:,2);
                        end
                    end

                    % score
                    roi = s.caTSA.roiArray{fi}.getRoiById(s.caTSA.ids(ci));
                    cci = find(dat.ids == s.caTSA.ids(ci));
                    if (length(cci) == 1 & ~isempty(roi))
                        dat.day(f).red_score(cci) = nanmean(dim{fi}(roi.indices));
                    end
                catch me
                    disp('REDNESS FAILED!');
                end
            end
        end
    end

% this is special case for the TWM which is a bit weird and complicated
function dat = get_two_whisker_model_fields (s, dat, f)
    fields_to_store = {'w1ExclusiveTouchTrials','w2ExclusiveTouchTrials','w1ProExclusiveTouchTrials','w1RetExclusiveTouchTrials', ...
                       'w2ProExclusiveTouchTrials','w2RetExclusiveTouchTrials'};
    m_fname = sprintf('dffBased_%sAbsMaxKappaZeroNotouch-%sAbsMaxKappaZeroNotouchtwoWhiskerModelWithCrossvalScore', s.whiskerTag{1}, s.whiskerTag{2});

    % for rise fall time
    m_datfit_fname  = sprintf('dffBased_%sAbsMaxKappaZeroNotouch-%sAbsMaxKappaZeroNotouchtwoWhiskerModelAllDataFit', s.whiskerTag{1}, s.whiskerTag{2});
    dat.day(f).w1_tRise = nan*dat.ids;
    dat.day(f).w1_tDecay = nan*dat.ids;
    dat.day(f).w2_tRise = nan*dat.ids;
    dat.day(f).w2_tDecay = nan*dat.ids;    

    if (isfield(s.caTSAArray,'caTSA'))
        for fi=1:length(fields_to_store)
            dat.day(f).(fields_to_store{fi}) = nan*dat.ids;
        end
        for ci=1:length(s.caTSAArray.caTSA)
            if (~isempty(s.caTSAArray.caTSA{ci}.cellFeatures.get(m_fname)))
                ti = find(ismember(s.caTSAArray.caTSA{ci}.ids, dat.ids));
                di = find( ismember(dat.ids,s.caTSAArray.caTSA{ci}.ids));

                m = struct_nan_blanks(s.caTSAArray.caTSA{ci}.cellFeatures.get(m_fname));
                for fi=1:length(fields_to_store)
                    fv = [m.(fields_to_store{fi})];
                    dat.day(f).(fields_to_store{fi})(di) = fv(ti);
                end 

                m_datfit = struct_nan_blanks(s.caTSAArray.caTSA{ci}.cellFeatures.get(m_datfit_fname));
                for i=1:length(s.caTSAArray.caTSA{ci}.ids)
                    ii = find(dat.ids == s.caTSAArray.caTSA{ci}.ids(i));
                    if (~isempty(ii) && isstruct(m_datfit) && isfield(m_datfit, 'w1') && isstruct(m_datfit(i).w1))
                        dat.day(f).w1_tRise(ii) = m_datfit(i).w1.tRise;
                        dat.day(f).w1_tDecay(ii) = m_datfit(i).w1.tDecay;
                        dat.day(f).w2_tRise(ii) = m_datfit(i).w2.tRise;
                        dat.day(f).w2_tDecay(ii) = m_datfit(i).w2.tDecay;
                    end
                end
            end
        end           
    else % caTSA only
        for fi=1:length(fields_to_store)
            dat.day(f).(fields_to_store{fi}) = nan*dat.ids;
        end
        if (~isempty(s.caTSA.cellFeatures.get(m_fname)))
            ti = find(ismember(s.caTSA.ids, dat.ids));
            di = find( ismember(dat.ids,s.caTSA.ids));

            m = struct_nan_blanks(s.caTSA.cellFeatures.get(m_fname));
            for fi=1:length(fields_to_store)
                fv = [m.(fields_to_store{fi})];
                dat.day(f).(fields_to_store{fi})(di) = fv(ti);
            end

            m_datfit = struct_nan_blanks(s.caTSA.cellFeatures.get(m_datfit_fname));
            for i=1:length(s.caTSA.ids)
                ii = find(dat.ids == s.caTSA.ids(i));
                if (~isempty(ii) && isstruct(m_datfit) && isfield(m_datfit, 'w1') && isstruct(m_datfit(i).w1))
                    dat.day(f).w1_tRise(ii) = m_datfit(i).w1.tRise;
                    dat.day(f).w1_tDecay(ii) = m_datfit(i).w1.tDecay;
                    dat.day(f).w2_tRise(ii) = m_datfit(i).w2.tRise;
                    dat.day(f).w2_tDecay(ii) = m_datfit(i).w2.tDecay;
                end
            end
        end
    end

% works for cases where we just have a vector with one value per ID
function dat = add_single_evoked_feature(s, dat, f, field_name, wh_replace) 
    if (nargin < 5 || isempty(wh_replace)) ; wh_replace = 0 ; end

    if (wh_replace)
        if (length(s.whiskerTag) >= 1)
            field_name_out = strrep(field_name, s.whiskerTag{1}, 'W1');;
        end
        if (length(s.whiskerTag) >= 2)
            field_name_out = strrep(field_name_out, s.whiskerTag{2}, 'W2');;
        end
    else
        field_name_out = field_name;
    end

    if (isfield(s.caTSAArray,'caTSA'))
        dat.day(f).(field_name_out) = nan*dat.ids;
        for ci=1:length(s.caTSAArray.caTSA)
            ti = find(ismember(s.caTSAArray.caTSA{ci}.ids, dat.ids));
            di = find( ismember(dat.ids,s.caTSAArray.caTSA{ci}.ids));
            if (~isempty(s.caTSAArray.caTSA{ci}.cellFeatures.get(field_name)))
                val = s.caTSAArray.caTSA{ci}.cellFeatures.get(field_name);
                dat.day(f).(field_name_out)(di) = val(ti);
            end
        end
    else
        dat.day(f).(field_name_out) = nan*dat.ids;
        ti = find(ismember(s.caTSA.ids, dat.ids));
        di = find( ismember(dat.ids,s.caTSA.ids));
        if (~isempty(s.caTSA.cellFeatures.get(field_name)))
            val = s.caTSA.cellFeatures.get(field_name);
            dat.day(f).(field_name_out)(di) = val(ti);
        end
    end

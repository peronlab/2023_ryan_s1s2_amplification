%
% Aggregator that will grab all s1s2 lesion data
%
function all_dat = get_s1s2_all_dat (force_redo)
    if (nargin < 1) ; force_redo = 0 ;end

    % thresholds for type determination
    ev_touch_thresh = 0.1; 
    twm_touch_thresh = 0.1; 
    ev_whisking_thresh = 0.05;
    glm_whisking_thresh = 0.2;

    ev_touch_thresh = 0.05; 
    ev_touch_thresh = 0.05; 
   

    % path stuff -- where is the data -- this is based on hostname
    [irr hname] = system('hostname');
    hname = hname(1:end-1);

    switch hname
        case 'apollo' % CLONE THIS SECTION with your own hostname
            root_dir = '/Volumes/Data/s1s2_lesion/';
            
        case 'Laurens-MacBook-Pro-2.local' %LR laptop
            %root_dir = '~/Desktop/poster_tmp/';%temp-- need new harddrive. cl
            root_dir= '/Volumes/volume_imaging/'

        otherwise
            disp(['HOSTNAME ' hname ' NOT RECOGNIZED; YOU WILL HAVE A PROBLEM']);
    end

    % load animal data - go directory by directory
    mdat = get_s1s2_metadata;
    mdat.root_dir = root_dir;
    all_dat = mdat;
    opwd = pwd;
    for a=1:length(mdat.anims)
        disp(['Processing ' mdat.anims{a}]);    

        % TEMPORARY - we only care about the volume mice for now
%   if (~mdat.has_deep_data_s1s2(a)) ; disp('FOR NOW skipping animals that are not 2 subvolumes per area') ; continue ; end

        % call to get_s1s2_single_animal_data
        cd([root_dir mdat.anims{a} '/session_neuropilone']); 
        all_dat.anim_data{a} = get_s1s2_single_animal_data(force_redo);
       
        % build matrices for desired variables -- fields start out with these names
        src_fn = {'evokedDffScores_probResponseTouchW1ExcludeW2','evokedDffScores_meanResponseTouchW1ExcludeW2', ...
                  'evokedDffScores_probResponseTouchW2ExcludeW1', 'evokedDffScores_meanResponseTouchW2ExcludeW1', ...
                  'evokedDffScores_probResponseTouchW1pExclusive','evokedDffScores_meanResponseTouchW1pExclusive', ...
                  'evokedDffScores_probResponseTouchW1rExclusive','evokedDffScores_meanResponseTouchW1rExclusive', ...
                  'evokedDffScores_probResponseTouchW2pExclusive','evokedDffScores_meanResponseTouchW2pExclusive', ...
                  'evokedDffScores_probResponseTouchW2rExclusive','evokedDffScores_meanResponseTouchW2rExclusive', ...
                  'evokedDffScores_probResponseWhiskExcludeOptogenTouch','evokedDffScores_meanResponseWhiskExcludeOptogenTouch', ...
                  'dffBased_k16t14_MeanWhiskerAmplitudeExcludeTouchGLMScore','dffBased_k16t14_WhiskerSetpointExcludeTouchGLMScore', ...
                  'w1ProExclusiveTouchTrials', 'w1RetExclusiveTouchTrials', 'w2ProExclusiveTouchTrials', 'w2RetExclusiveTouchTrials'}; 
                  %'dffBased_k16t14_MeanWhiskerAmplitudeGLMScore','dffBased_k16t14_WhiskerSetpointGLMScore', ...

        % we name the fields this to make names shorter and easier
        fn = {'probRespW1', 'meanRespW1', ... % all of these are 'evoked' features
              'probRespW2', 'meanRespW2', ...
              'probRespW1p', 'meanRespW1p', ...
              'probRespW1r', 'meanRespW1r', ...
              'probRespW2p', 'meanRespW2p', ...
              'probRespW2r', 'meanRespW2r', ...
              'probRespWhisk', 'meanRespWhisk', ...
              'whAmplitude', 'whSetpoint', ...
              'w1p', 'w1r', 'w2p', 'w2r'}; % two-whisker model output

        % build blank matrices that are nDay X nNeurons long
        for f=1:length(fn)
            all_dat.anim_mats{a}.(fn{f}) = nan*zeros(length(all_dat.anim_data{a}.day), length(all_dat.anim_data{a}.ids));
        end

        % some stuff to make life easier
        date_str = {};
        is_merged_across_days = zeros(1,length(all_dat.anim_data{a}.day));
        for d=1:length(all_dat.anim_data{a}.day) ; 
            sfn = all_dat.anim_data{a}.day(d).source_filename;
            ui = find(sfn == '_');
            if (length(ui) == 4) % this means we have a real date; pull MM_DD
                date_str{d} = sfn(ui(2)+1:ui(4)-1); % MM_DD
            elseif (length(ui) == 3)
                date_str{d} = sfn(ui(1)+1:ui(3)-1); % vol_???
                is_merged_across_days(d) = 1; 
            else
                date_str{d} = ''; % DEGENERATE!
            end
        end 
        all_dat.anim_data{a}.date_str = date_str;
        all_dat.anim_data{a}.is_merged_across_days = is_merged_across_days;

        % loop over all sessions and get data
        has_data = ones(1,length(all_dat.anim_data{a}.day));
        for d=1:length(all_dat.anim_data{a}.day);
            if (isempty(all_dat.anim_data{a}.day(d).source_filename))
                disp(sprintf('DATA EMPTY: (a,d) (%d, %d), which usually is OK', a, d));
                has_data(d) = 0;
            end
            for f=1:length(fn)
                vv = nan*all_dat.anim_data{a}.ids;
                if (~isempty(all_dat.anim_data{a}.day(d).(src_fn{f})) && length(find(isnan(all_dat.anim_data{a}.day(d).(src_fn{f})))) ~= length(all_dat.anim_data{a}.day(d).(src_fn{f})))
                    vv = all_dat.anim_data{a}.day(d).(src_fn{f});
                elseif (has_data(d))
                    disp(sprintf('DATA MISSING: (a,d) (%d, %d): %s field: %s', a, d, all_dat.anim_data{a}.day(d).source_filename, src_fn{f}));
                end
                all_dat.anim_mats{a}.(fn{f})(d,:) = vv;
            end    
        end
        all_dat.anim_data{a}.has_data = has_data;

        % use threshold to assign types -- we have this so we only have ONE place where we decide if a cell is or is not touch/whisking/etc.
        for d=1:length(all_dat.anim_data{a}.day)
            % most basic 
            all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d} = find(all_dat.anim_mats{a}.probRespW1p(d,:) > ev_touch_thresh);
            all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d} = find(all_dat.anim_mats{a}.probRespW1r(d,:) > ev_touch_thresh);
            all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d} = find(all_dat.anim_mats{a}.probRespW2p(d,:) > ev_touch_thresh);
            all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d} = find(all_dat.anim_mats{a}.probRespW2r(d,:) > ev_touch_thresh);
if (0)
            ev_touch_thresh = 0.1;
            all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d} = find(all_dat.anim_mats{a}.meanRespW1p(d,:) > ev_touch_thresh);
            all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d} = find(all_dat.anim_mats{a}.meanRespW1r(d,:) > ev_touch_thresh);
            all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d} = find(all_dat.anim_mats{a}.meanRespW2p(d,:) > ev_touch_thresh);
            all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d} = find(all_dat.anim_mats{a}.meanRespW2r(d,:) > ev_touch_thresh);
end
if (0) % display numbbers so you can make sure threshold is not yielding ridiculous low/high neuron count
adt =  all_dat.types_by_idx{a};
nc = length(all_dat.anim_data{a}.ids);
disp([all_dat.anims{a} ' day: ' num2str(d) ' fraction : ' num2str([length(adt.ev_touch_w1p_by_day{d}) length(adt.ev_touch_w1r_by_day{d}) length(adt.ev_touch_w2p_by_day{d}) length(adt.ev_touch_w2r_by_day{d}) ]/nc)]);
disp([all_dat.anims{a} ' day: ' num2str(d) ' n: ' num2str([length(adt.ev_touch_w1p_by_day{d}) length(adt.ev_touch_w1r_by_day{d}) length(adt.ev_touch_w2p_by_day{d}) length(adt.ev_touch_w2r_by_day{d}) nc])]);
end
            all_dat.types_by_idx{a}.twm_touch_w1p_by_day{d} = find(all_dat.anim_mats{a}.w1p(d,:) > twm_touch_thresh);
            all_dat.types_by_idx{a}.twm_touch_w1r_by_day{d} = find(all_dat.anim_mats{a}.w1r(d,:) > twm_touch_thresh);
            all_dat.types_by_idx{a}.twm_touch_w2p_by_day{d} = find(all_dat.anim_mats{a}.w2p(d,:) > twm_touch_thresh);
            all_dat.types_by_idx{a}.twm_touch_w2r_by_day{d} = find(all_dat.anim_mats{a}.w2r(d,:) > twm_touch_thresh);

            all_dat.types_by_idx{a}.ev_whisking_by_day{d} = find(all_dat.anim_mats{a}.probRespWhisk(d,:) > ev_whisking_thresh);
            all_dat.types_by_idx{a}.glm_whisking_amp_by_day{d} = find(all_dat.anim_mats{a}.whAmplitude(d,:) > glm_whisking_thresh);
            all_dat.types_by_idx{a}.glm_whisking_sp_by_day{d} = find(all_dat.anim_mats{a}.whSetpoint(d,:) > glm_whisking_thresh);

            % bidirectional / sw / mw, based on either evoked or two whisker model (if statement below)
            touch_mat = zeros(4,size(all_dat.anim_mats{a}.probRespW1,2));
            if (1) % use event-based
                touch_mat(1,all_dat.types_by_idx{a}.ev_touch_w1p_by_day{d}) = 1;
                touch_mat(2,all_dat.types_by_idx{a}.ev_touch_w1r_by_day{d}) = 1;
                touch_mat(3,all_dat.types_by_idx{a}.ev_touch_w2p_by_day{d}) = 1;
                touch_mat(4,all_dat.types_by_idx{a}.ev_touch_w2r_by_day{d}) = 1;
            else % use model
                touch_mat(1,all_dat.types_by_idx{a}.twm_touch_w1p_by_day{d}) = 1;
                touch_mat(2,all_dat.types_by_idx{a}.twm_touch_w1r_by_day{d}) = 1;
                touch_mat(3,all_dat.types_by_idx{a}.twm_touch_w2p_by_day{d}) = 1;
                touch_mat(4,all_dat.types_by_idx{a}.twm_touch_w2r_by_day{d}) = 1;
            end
            all_dat.types_by_idx{a}.usw_by_day{d} = find(sum(touch_mat) == 1);
            all_dat.types_by_idx{a}.bsw_by_day{d} = union(find(sum(touch_mat(1:2,:)) == 2 & sum(touch_mat(3:4,:)) == 0), ...
                                                          find(sum(touch_mat(1:2,:)) == 0 & sum(touch_mat(3:4,:)) == 2));
            w1i = find(sum(touch_mat(1:2,:)) > 0); 
            w2i = find(sum(touch_mat(3:4,:)) > 0); 
            all_dat.types_by_idx{a}.mw_by_day{d} = intersect(w2i,w1i);

            all_dat.types_by_idx{a}.touch_w1p_by_day{d} = find(touch_mat(1,:) == 1);
            all_dat.types_by_idx{a}.touch_w1r_by_day{d} = find(touch_mat(2,:) == 1);
            all_dat.types_by_idx{a}.touch_w2p_by_day{d} = find(touch_mat(3,:) == 1);
            all_dat.types_by_idx{a}.touch_w2r_by_day{d} = find(touch_mat(4,:) == 1);
        end

        % redness ... 
        if (mdat.retrograde_label_type(a) > 0)
            red_score = all_dat.anim_data{a}.day(mdat.redness_day_idx(a)).red_score;
            all_dat.types_by_idx{a}.is_red = find(red_score >mdat.redness_thresh(a));
        end

        % xyz
        invali = find(isnan(all_dat.anim_data{a}.z_um));
        s1i = find(all_dat.anim_data{a}.ids >= all_dat.s1_id_range{a}(1) & all_dat.anim_data{a}.ids <= all_dat.s1_id_range{a}(2));
        s1_sv = unique(floor(all_dat.anim_data{a}.ids(s1i)/10000000));
        s2i = find(all_dat.anim_data{a}.ids >= all_dat.s2_id_range{a}(1) & all_dat.anim_data{a}.ids <= all_dat.s2_id_range{a}(2));
        s2_sv = unique(floor(all_dat.anim_data{a}.ids(s2i)/10000000));
        plni = rem(floor(all_dat.anim_data{a}.ids(:)/10000),10);
        disp([mdat.anims{a} ': ' num2str(length(invali)) ' invalid cells; s1: ' num2str(length(s1i)) ' s2: ' num2str(length(s2i))]);
        for i=1:length(invali)
            y_frac = all_dat.anim_data{a}.y_um(invali(i))/600;
            if (ismember(invali(i), s1i))
                base_z = mdat.s1_start_depth_super_rel_L1_micron(a);
                svi = find(s1_sv == floor(all_dat.anim_data{a}.ids(invali(i))/10000000));
                if (svi > 1) ; base_z = base_z+(svi-1)*3*mdat.dz_used(a); end
            else
                base_z = mdat.s2_start_depth_super_rel_L1_micron(a);
                svi = find(s2_sv == floor(all_dat.anim_data{a}.ids(invali(i))/10000000));
                if (svi > 1) ; base_z = base_z+(svi-1)*3*mdat.dz_used(a); end                
            end
            all_dat.anim_data{a}.z_um(invali(i)) = base_z + y_frac*mdat.dz_used(a) + (plni(invali(i))-1)*mdat.dz_used(a);
        end

    end

    cd (opwd);


function idx = get_top_indices(vv, q_thresh)
    vv (find(isnan((vv))))= -1;
    [irr sorti] = sort(vv, 'descend');
    idx = sorti(1:ceil(q_thresh*length(sorti))); 


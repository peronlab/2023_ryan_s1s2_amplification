function stats = get_session_beh_stats (s)

    %% basic stuff on correct/incorrect/touched etc.

    n_trials = length(s.trialIds);

    % lick info for touch restriction to prelick
    lickESL = s.behavESA.getEventSeriesByIdStr('Left lick',1).copy();
    lickESR = s.behavESA.getEventSeriesByIdStr('Right lick',1).copy();
    allLickTimes = union(lickESL.eventTimes, lickESR.eventTimes);
    allLickES = pldo.eventSeries();
    allLickES.id = 1; 
    allLickES.idStr = 'All Licks';
    allLickES.type = 1;
    allLickES.eventTimes = allLickTimes;
    allLickES.updateEventSeriesFromTrialTimes(s.behavESA.trialTimes);

    % touch trials - ONLY IF TOUCH IS BEFORE FIRST LICK
    t_trials = unique(s.whiskerBarContactESA.esa{1}.eventTrials);
    bad_t_trials = [];
    for t=1:length(t_trials)
        touch_t = s.whiskerBarContactESA.esa{1}.eventTimes(find(s.whiskerBarContactESA.esa{1}.eventTrials == t_trials(t)));
        lti = find(allLickES.eventTrials == t_trials(t));
        if (length(lti) > 0)
            first_lick_t = allLickES.eventTimes(lti(1));
            touch_t = touch_t(find(touch_t < first_lick_t));
            if (isempty(touch_t)) ; bad_t_trials = [bad_t_trials t_trials(t) ] ; end
        end
    end
    t_trials = setdiff(t_trials, bad_t_trials);
    touchi = find(ismember(s.trialIds, t_trials));
    non_touchi = find(~ismember(s.trialIds, t_trials));

   
    % touch epochs isolate - ALL touches and filtered below to pre-lick
    baseTime = s.whiskerAngleTSA.time;
    dt = mode(diff(baseTime));

    proContactTS =  s.whiskerBarContactClassifiedESA.esa{1}.deriveTimeSeries(baseTime, s.whiskerAngleTSA.timeUnit, [0 1]); % gives you a binary vector of length baseTime with 1 where touches occurred, 0 not
    retContactTS =  s.whiskerBarContactClassifiedESA.esa{2}.deriveTimeSeries(baseTime, s.whiskerAngleTSA.timeUnit, [0 1]);
    proTouchIdx = find(proContactTS.value);
    retTouchIdx = find(retContactTS.value);
    touchIdx = union(proTouchIdx, retTouchIdx);

    uProTrials = unique(s.whiskerCurvatureChangeTSA.trialIndices(1,proTouchIdx));
    netProKappa = nan*zeros(1,length(uProTrials));
    for t=1:length(uProTrials)
        ii = find(s.whiskerCurvatureChangeTSA.trialIndices == uProTrials(t));
        ii = intersect(ii,proTouchIdx);

        netProKappa(t) = nansum(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,ii)));
    end

    uRetTrials = unique(s.whiskerCurvatureChangeTSA.trialIndices(1,proTouchIdx));
    netRetKappa = nan*zeros(1,length(uRetTrials));
    for t=1:length(uRetTrials)
        ii = find(s.whiskerCurvatureChangeTSA.trialIndices == uRetTrials(t));
        ii = intersect(ii,retTouchIdx);

        netRetKappa(t) = nansum(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,ii)));
    end    


    %% Hilbert transform - look to Hill and Kleinfeld 2011, Neuron, Figure 3 for a detailed exposition
    % max theta, velocity (99.9th %ile - max is degenerate many times)

    % hilbert ; max amplitude, mean setpoint 
    whAngTS = s.whiskerAngleTSA.getTimeSeriesByIdx(1);
    whDtSec = pldo.timeSeries.convertTime(mode(diff(whAngTS.time)), whAngTS.timeUnit, pldo.timeSeries.second);
    whAngVec = whAngTS.value;
    whAngVec(find(isnan(whAngVec))) = 0; % otherwise this breaks

    sampleRate=  1/whDtSec;
    BandPassCutOffsInHz = [6 30];  %%check filter parameters!!!
    W1 = BandPassCutOffsInHz(1) / (sampleRate/2);
    W2 = BandPassCutOffsInHz(2) / (sampleRate/2);
    [b,a]=butter(2,[W1 W2]);
    filteredSignal = filtfilt(b, a, whAngVec);

    [b,a]=butter(2, 6/ (sampleRate/2),'low');
    setpoint = filtfilt(b,a,whAngVec-filteredSignal);

    hh=hilbert(filteredSignal);
    amplitude = abs(hh);

    velocity = diff(whAngVec)/whDtSec;
    velocity(end+1) = nan; % match length

    % per trial stuff -- pole position, etc.
    trial_dat = [];


    s.generateContactPropertiesHashes;
    theta_at_touch = s.whiskerBarContactClassifiedESA.esa{1}.eventPropertiesHash.get('thetaAtTouchOnset');
    for t=1:length(s.trial)
        if (ismember(s.trial{t}.typeIds, [1 2]))
            trial_dat(t).correct = 1; %1: correct; 0: incorrect ; -1: ignore
        elseif (ismember(s.trial{t}.typeIds, [3 4]))
            trial_dat(t).correct = 0; %1: correct; 0: incorrect ; -1: ignore
        else
            trial_dat(t).correct = -1; %1: correct; 0: incorrect ; -1: ignore
        end

        if (ismember(s.trial{t}.typeIds, [2 3]))
            trial_dat(t).lick_right = 1; %1: right; 0: left ; -1: ignore
        elseif (ismember(s.trial{t}.typeIds, [1 4]))
            trial_dat(t).lick_right = 0; %1: right; 0: left ; -1: ignore
        else
            trial_dat(t).lick_right = -1; %1: right; 0: left ; -1: ignore
        end

        trial_dat(t).pole_position = s.trial{t}.behavParams.get('polePositionAP');

        % COMPUTE NET dK, and touch count 
        %   - mean amplitude
        %   - max amplitude
        %   - mean setpoint
        %   - max setpoint
        %   - max velocity
        trial_dat(t).net_kappa = nan;
        trial_dat(t).net_kappa_pre_first_lick = nan;
        trial_dat(t).q99_kappa_pre_first_lick = nan;
        trial_dat(t).num_pro_touches = 0;
        trial_dat(t).num_pro_touches_pre_first_lick = 0;
        trial_dat(t).num_ret_touches = 0;
        trial_dat(t).num_ret_touches_pre_first_lick = 0;
        trial_dat(t).mean_theta_at_touch = nan; % ONLY FOR PROTRACTIONS
        trial_dat(t).angle_distro = [];

        trial_dat(t).mean_amplitude_pre_first_lick = nan;
        trial_dat(t).q99_amplitude_pre_first_lick = nan;
        trial_dat(t).mean_setpoint_pre_first_lick = nan;
        trial_dat(t).q99_setpoint_pre_first_lick = nan;
        trial_dat(t).mean_velocity_pre_first_lick = nan;
        trial_dat(t).q99_velocity_pre_first_lick = nan;
    
        lick_this_triali = find(allLickES.eventTrials == s.trial{t}.id);
        if (length(lick_this_triali) > 0)
            t_first_lick = min(allLickES.eventTimes(lick_this_triali));
        else
            t_first_lick = Inf;
        end
    
        % kappa based stuff
        if (ismember(s.trial{t}.id, uProTrials))
            pti = intersect(find(s.whiskerCurvatureChangeTSA.trialIndices == s.trial{t}.id), proTouchIdx);
            trial_dat(t).net_kappa = nansum(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,pti)));

            pti_pre_lick = pti(find(s.whiskerCurvatureChangeTSA.time(pti) < t_first_lick));
            trial_dat(t).net_kappa_pre_first_lick = nansum(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,pti_pre_lick)));
            if (length(pti_pre_lick) > 0)
                trial_dat(t).q99_kappa_pre_first_lick = quantile(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,pti_pre_lick)), .99);
            end

            pro_touchi = find(s.whiskerBarContactClassifiedESA.esa{1}.eventTrials == s.trial{t}.id);
            if (length(pro_touchi) > 0)
                pro_touchi = pro_touchi(1:2:end); % only consider touch *start* time
                trial_dat(t).num_pro_touches = length(pro_touchi);
                trial_dat(t).num_pro_touches_pre_first_lick = length(find(s.whiskerBarContactClassifiedESA.esa{1}.eventTimes(pro_touchi) < t_first_lick));

                % collect touch angles
                trial_dat(t).mean_theta_at_touch = nanmean(theta_at_touch(pro_touchi));
            end
        end
        if (ismember(s.trial{t}.id, uRetTrials))
            rti = intersect(find(s.whiskerCurvatureChangeTSA.trialIndices == s.trial{t}.id), retTouchIdx);
            trial_dat(t).net_kappa = trial_dat(t).net_kappa+nansum(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,rti)));

            rti_pre_lick = rti(find(s.whiskerCurvatureChangeTSA.time(rti) < t_first_lick));
            trial_dat(t).net_kappa_pre_first_lick = trial_dat(t).net_kappa_pre_first_lick + nansum(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,rti_pre_lick)));
            if (length(rti_pre_lick) > 0)
                trial_dat(t).q99_kappa_pre_first_lick = trial_dat(t).q99_kappa_pre_first_lick + quantile(abs(s.whiskerCurvatureChangeTSA.valueMatrix(1,rti_pre_lick)), .99);
            end

            ret_touchi = find(s.whiskerBarContactClassifiedESA.esa{2}.eventTrials == s.trial{t}.id);
            if (length(ret_touchi) > 0)
                ret_touchi = ret_touchi(1:2:end); % only consider touch *start* time
                trial_dat(t).num_ret_touches = length(ret_touchi);
                trial_dat(t).num_ret_touches_pre_first_lick = length(find(s.whiskerBarContactClassifiedESA.esa{2}.eventTimes(ret_touchi) < t_first_lick));
            end
        end

        % get whisker angle values from pole in reach up to either first touch, first lick, or pole out of reach
        bti = find(s.whiskerBarInReachES.eventTrials == s.trial{t}.id);
        if (length(bti) == 2)
            start_ti = min(find(s.whiskerAngleTSA.time > s.whiskerBarInReachES.eventTimes(bti(1))));
            touchi = find(s.whiskerBarContactESA.esa{1}.eventTrials == s.trial{t}.id);
            if (length(touchi) > 0)
                t_first_touch = min(s.whiskerBarContactESA.esa{1}.eventTimes(touchi));
            else
                t_first_touch = Inf;
            end

            end_time = min(t_first_lick, s.whiskerBarInReachES.eventTimes(bti(2)));
            end_time = min(end_time, t_first_touch);
            end_ti = max(find(s.whiskerAngleTSA.time < end_time));

            if (length(start_ti) == 1 & length(end_ti) == 1)
                valti = find(s.whiskerAngleTSA.time <= end_ti & s.whiskerAngleTSA.time >= start_ti);

                if (length(valti) > 0)
                    trial_dat(t).angle_distro = s.whiskerAngleTSA.valueMatrix(1,valti);

                    trial_dat(t).mean_amplitude_pre_first_lick = nanmean(amplitude(valti));
                    trial_dat(t).q99_amplitude_pre_first_lick = quantile(amplitude(valti), .99);
                    trial_dat(t).mean_setpoint_pre_first_lick = nanmean(amplitude(valti));
                    trial_dat(t).q99_setpoint_pre_first_lick = quantile(setpoint(valti), .99);
                    trial_dat(t).mean_velocity_pre_first_lick = nanmean(velocity(valti));
                    trial_dat(t).q99_velocity_pre_first_lick = quantile(velocity(valti), .99);
                end
            end
        end    

        % other whisker stuff: for pole in reach period, terminating on first lick if earlier
        %   - mean amplitude
        %   - max amplitude
        %   - mean setpoint
        %   - max setpoint
        
    end
    stats.trial_dat = trial_dat;
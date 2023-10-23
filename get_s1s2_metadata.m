function mdat = get_s1s2_metadata
    % --- actual metadata
    mdat.anims = {'an014351', 'an017518', 'an014359', 'an017722', 'an018921', 'an014332', 'an018490', 'an018489', 'an018920', 'an020047', 'an020428',  ... %vs1 lesion
        'an016623', 'an014362', 'an017517', 'an017721', 'an018927', 'an018933', 'an020048', 'an018922',... %vs2 lesion
        'an016650', 'an016652','an014363','an017510', ... %s1s2 only
        'an014353', 'an014355', 'an014354', 'an014335', ... %1d sham only (some have CM too, but deal w.later)
        'an018912', 'an018919', 'an020055', 'an019781', 'an019770','an020035', 'an020044', 'an020421', 'an019806', 'an020430', 'an019776', 'an019821'}; % retro label only (though some have sham!)
    
    mdat.valid_dates={{'06_26','06_27','07_01','07_02'},{'08_22','08_23','08_24','08_26'},{'08_25', '08_29'},{'09_18','09_19','09_20','09_21'}, {'10_18','10_19','10_20','10_21'}, {'10_27','10_31','11_01'},{'10_21','10_25','10_27','10_28'}, {'11_18','11_22','11_23','11_24'}, {'11_18','11_22','11_23','11_24'}, {'02_06'},{'02_06'},  ...
        {'07_29','08_02','08_05','08_10','08_11'}, {'08_01','08_05','08_10'},{'08_22', '08_23','08_26','08_29'},{'10_03','10_04','10_05','10_06','10_11','10_12'},{'10_11','10_12','10_13','10_14','10_15','10_16'},{'10_14','10_15','10_16','10_17','10_18','10_19'},{'11_25', '11_28','11_29','11_30'},{'11_25', '11_28','11_29','11_30'}, ...
        {'05_26','05_27','05_31'},{'06_03','06_08'},{'08_08','08_09'},{'08_17','08_19'} ...
        {'09_18','09_19','09_20','09_21'},{'09_11','09_12','09_13', '09_14'},{'09_11','09_12','09_13', '09_14'}, {'10_20','10_24','10_25','10_26'}, ...
        {'12_13'},{'12_13'},{'12_14'},{'01_15','01_18'},{'01_15','01_18'},{'01_17'},{'02_07'},{'02_13'},{'02_15'},{'02_15'},{'02_24'},{'02_27'}};

    % -1 = use for pre / 1 = use for post sham ; usually use 2 pre days but sometimes day 1 is weird (1st day on img rig) ; 1 post day
    mdat.sham_dates={[0 0 0 0],[0 0 0 0],[0 0], [0 0 0 0], [0 0 0 0], [0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0], [0], [0],  ...
        [-1 -1 1 0 0], [-1 1 0],[0 -1 1 0], [-1 -1 1 0 0 0 ], [-1 -1 1 0 0 0], [-1 -1 1 1 0 0], [0 0 0 0],[0 0 0 0],  ...
        [ 0 0 0],[ 0 0],[0 0], [ 0 0], ...
        [0 -1 1 0],[-1 0 1 0], [-1 -1 1 1], [0 -1 1 0], ...
        [0], [0], [0], [-1 1], [-1 1],[0], [0], [0], [0],[0],[0],[0]};

    % -1 = use for pre / 1 = use for post lesion ; usually use 2 pre days but sometimes day 1 is weird (1st day on img rig) ; 1 post day
    mdat.lesion_dates={[-1 -1 1 0],[0 -1 -1 1],[-1 1], [-1 -1 1 1], [-1 -1 1 0 ], [-1 -1 1], [-1 -1 1 0], [-1 0 1 0], [-1 0 1 0] ,  [0], [0],  ...
        [0 0 -1 1 0], [0 -1 1],[0 0 -1 1], [0 0 0 0 -1 1], [0 0 -1 -1 1 0], [0 0 -1 1 0 ], [-1 -1 1 0], [-1 -1 1 0], ...
        [ 0 0 0],[ 0 0],[0 0], [ 0 0], ...
        [0 0 0 0],[-1 -1 1 1], [0 0 0 0], [0 0 0 0], ...
        [0], [0], [0], [0 0], [0 0], [0], [0], [0], [0],[0],[0],[0]};

    %% DAY 2: UNCOMMENT
%    mdat.lesion_dates={[-1 -1 0 1],[0 -1 -1 1],[-1 1], [-1 -1 0 1], [-1 -1 0 1 ], [-1 -1 1], [-1 -1 0 1], [-1 0 0 1], [-1 0 0 1] ,  [0], [0],  ...
%        [0 0 -1 0 1], [0 -1 1],[0 0 -1 1], [0 0 0 0 -1 1], [0 0 -1 -1 0 1], [0 0 -1 0 1 ], [-1 -1 0 1], [-1 -1 0 1], ...
%        [ 0 0 0],[ 0 0],[0 0], [ 0 0], ...
%        [0 0 0 0],[-1 -1 0 1], [0 0 0 0], [0 0 0 0], ...
%        [0], [0], [0], [0 0], [0 0], [0], [0], [0], [0],[0],[0],[0]};

    % 1 = use for retrograde analysis (pre ANY lesion)
    mdat.retro_dates={[0 0 0 0],[0 0 1 0],[1 0], [0 0 0 0], [0 0 0 0], [0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0], [1], [1],   ...
        [0 0 0 0 0], [1 0 0],[0 0 0 0], [0 0 0 0 0 0 ], [0 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0], [0 0 0], ...
        [ 0 0 0],[ 0 0],[0 0], [ 0 0], ...
        [0 0 0 0],[1 0 0 0], [1 0 0 0], [0 0 0 0], ...
        [1], [1], [1], [1 0], [1 0], [1], [1], [1], [1],[1],[1],[1]};        

    % 1 = use for volumetric pre-anything analysis (pre ANY lesion)
    mdat.static_dates={[1 1 0 0],[1 1 1 0],[1 0], [0 0 0 0], [0 0 0 0], [0 0 0], [0 0 0 0], [0 0 0 0], [0 0 0 0],  [0], [0],  ...
        [1 1 0 0 0], [0 0 0],[1 1 0 0], [0 0 0 0 0 0 ], [0 0 0 0 0 0], [0 0 0 0 0 0], [0 0 0], [0 0 0], ...
        [ 1 1 1],[ 1 1],[1 1], [ 1 1], ...
        [0 0 0 0],[0 0 0 0], [1 0 0 0], [0 0 0 0], ...
        [0], [0], [0], [0 0], [0 0], [0], [0], [0], [0],[0],[0],[0]};        
       
       
    mdat.s1_id_range = {[10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000],[10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000],...
        [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000],  [10000000 30000000], [10000000 30000000], ...
        [10000000 30000000], [900000000 930000000], [10000000 30000000], [10000000 30000000], ...
        [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], ...
        [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [0 1], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000], [10000000 30000000]}; 
        %[a b] specifying range for s1; can be much bigger but don't make smaller

    mdat.s2_id_range = {[1000000000 1030000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000],[900000000 930000000], [900000000 930000000],[900000000 930000000], [900000000 930000000],[900000000 930000000],  ...
        [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], ...
        [1000000000 1030000000], [9000000000 9030000000], [900000000 930000000], [900000000 930000000], ...
        [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], ...
        [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [10000000 30000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000],[900000000 930000000], [900000000 930000000]}; 
        % note: for animals who don't have s2 imaging or something, just specified the range anyways, jsut copied from another mouse
  
    % lesion analysis in a few animals excludes some s2 planes
    mdat.s1_lesion_id_range = mdat.s1_id_range;
    mdat.s2_lesion_id_range = {[1000000000 1030000000], [900000000 910039999], [900000000 910039999], [900000000 930000000], [900000000 930000000], [900000000 930000000],[900000000 930000000], [900000000 930000000],[900000000 930000000],  ...
        [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000],[900000000 930000000], [900000000 930000000],  ...
        [1000000000 1030000000], [9000000000 9030000000], [900000000 930000000], [900000000 930000000], ...
        [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], ...
        [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [10000000 30000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000], [900000000 930000000]}; 
        % note: for animals who don't have s2 imaging or something, just specified the range anyways, jsut copied from another mouse


    mdat.s1_start_depth_super_rel_L1_micron = [52 44 28 0 0 0 0 0 0 0 0 ...
    62 0 56 0 0 0 0 0 ...
    0 108 50 52 ...
    0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0]; % depth rel L1 for TOP of first subvolume for s1

    mdat.s1_start_depth_deep_rel_L1_micron =[224 224 200 0 0 0 0 0 0 0 0 ...
    228 0 236 0 0 0 0 0 ...
    0 244 236 226 ...
    0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0];  % depth rel L1 for TOP of deeper subvolume for s1

    mdat.s2_start_depth_super_rel_L1_micron = [40 40 60 0 0 0 0 0 0 0 0 ...
    92 0 14 0 0 0 0 0 ...
    0 94 40 32 ...
    0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0]; %depth rel L1 for TOP of first subvolume for s2

    mdat.s2_start_depth_deep_rel_L1_micron = [220 220 240 0 0 0 0 0 0 0 0 ...
    272 0 188 0 0 0 0 0 ...
    0 234 220 212 ...
    0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0]; % depth rel L1 for TOP of deeper subvolume for s1depth for TOP of first subvolume for s2


    mdat.dz_used= [60 60 60 60 60 60 60 60 60 30 40 ...
    60 60 60 60 60 60 60 60 ...
    60 60 60 60 ...
    60 60 60 60 ...
    30 40 30 40 30 30 40 30 30 30 40 30]; % dz for volume imaging
    
    mdat.needs_subselect= [0 1 1 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 ...
    0 0 0 0 ...
    0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0]; %might not use, but keeping track of which animals were noted as needing further subselection of planes or cells due to issues. 1=needs help, 0=fine.
    
    mdat.sham_lesion_date = {'','', '','', '', '', '', '', '', '', '', ...
        '08_03', '','08_24','10_04','10_12','', '', '', ...
        '', '', '', '', ...
        '09_19','','09_12', '10_24', ...
        '', '', '', '01_17', '01_17', '', '','', '', '', '', ''}; % MM_DD format
        
    mdat.is_1daysham= [0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 1 1 1 0 0 ...
    0 0 0 0 ...
    1 0 0  1 ... % 14355 and 14354 have very low T cell counts -> exclude
    0 0 0 1 1 0 0 0 0 0 0]; %0= no sham or sham thats 2+ days, so we can exclude later if we want. 1= 1 day sham
   
%    mdat.s1_lesion_date = {'06_30','08_24','08_26','09_19', '10_19', '10_31', '10_25', '11_22','11_22', ...
    mdat.s1_lesion_date = {'06_30','','08_26','09_19', '10_19', '10_31', '10_25', '11_22','11_22', '', '', ...
        '', '','','','','', '', '', ...
        '','','','' ...
        '','','','' ...
        '','','', '', '', '', '', '', '', '', '', ''}; % 18490 10_25 removed -- very odd vS2 (but vS1 ok) ; 17518 very few T cells weird

    mdat.s1_lesion_whisker = [3 1 2 2 2 2 2 2 2 0 0 ...
         0 0 0 0 0 0 0 0 ...
         0 0 0 0 ...
         0 0 0 0 ...
         0 0 0 0 0 0 0 0 0 0 0 0];

%        '08_08', '08_08','08_26','10_11','10_14','10_17', '11_28', '11_28', ...
    mdat.s2_lesion_date = {'','','','', '', '', '', '', '', '', '',  ...
        '08_08', '','','10_11','10_14','10_17', '11_28', '11_28', ...
        '','','','' ...
        '','','','' ...
        '','', '', '', '', '', '', '', '', '', '', ''};  % 14362(2nd) was not cellres targeted; exclude?
    
    mdat.retrograde_label_type = [0 0 1 0 0 0 0 0 0 1 1 ...
    0 2 0 0 0 0 0 0 ...
    0 0 0 0 ...
    0 1 1 0 ... % 14355 and 14354
    2 2 2 2 2 1 1 2 1 2 2 1];  % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
    %2 2 2 2 2 1 1 2 1 2 2 1];  % 0 = none 1 = retrograde injected into s1 2 = retrograde injected into s2
 
    mdat.redness_thresh = [nan 0.7 0.75 nan nan nan nan nan nan .75 .75 ...
    nan 0.725 nan nan nan nan nan nan ...
    nan nan nan nan ...
    nan 0.75 0.75 nan ...
    0.65 0.7 0.75 0.75 0.75 0.75 0.70 0.70 0.75 0.75 0.75 0.75];  % above this is red 
 
    mdat.redness_day_idx = [nan 3 1 nan nan nan nan nan nan 1 1 ...
    nan 1 nan nan nan nan nan nan  ...
    nan nan nan nan ...
    nan 1 1 nan ...
    1 1 1 1 1 1 1 1 1 1 1 1];  % which day to use for redness (in all_dat.anim_dat)
 
    mdat.has_deep_data_s1s2 = [1 1 1 0 0 0 0 0 0 0 0 ...
    1 0 1 0 0 0 0 0 ...
    1 1 1 1 ...
    0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0]; % 1 if there are two subvolumes for deep as well as superficial (if only superficial, 0) for both s1s2
    
    mdat.has_effect_radius_data = [0 1 1 1 1 1 1 1 1 0 0 ...
    0 0 0 0 0 0 0 0 ...
    0 0 0 0 ...
    0 1 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0]; 

    
    % "derived" metadata

    % --- settings
    
    % should we load merged?
    mdat.load_merged = 1;

    % colors
    mdat.w1_touch_color = [1 0 0];
    mdat.w2_touch_color = [0 0.5 1];
    mdat.whisking_color = [0 .8 0];

    mdat.mw_color = [255 128 255]/255;
    mdat.usw_color = [64 192 192]/255;
    mdat.bsw_color = [64 128 192]/255;
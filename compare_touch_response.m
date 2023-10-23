% this will compare touch response for two session objects
function compare_touch_response(s_pre, s_post, subvol_id)
    if (nargin < 3) ; subvol_id = []; end

    ani_str = fileparts_crossplatform(fileparts_crossplatform(pwd),2);

    % string?
    if (~isobject(s_pre))
        load(s_pre);
        s_pre = s;
        load(s_post);
        s_post = s;
    end

    % Which fields to use? can change to mean for example
    w1_field = 'evokedDffScores_probResponseTouchW1ExcludeW2';
    w2_field = 'evokedDffScores_probResponseTouchW2ExcludeW1';
    score_thresh = 0.15; % must be this to be considered touch cell

    % get touch scores 
    if (isfield(s_pre.caTSAArray,'caTSA')) % are we dealing with mutliple subvolumes?
        w1_t_pre = [];
        w2_t_pre = [];
        for c=1:length(s_pre.caTSAArray.caTSA) % this assumes pre/post are same size
            if (c == 1 && isempty( s_pre.caTSAArray.caTSA{c}.cellFeatures.get(w1_field))) ; s_pre.computeEvokedFeatures; end

            % matched subvolume?
            if (~isempty(subvol_id) & ~strcmp(subvol_id, s_pre.caTSAArray.caTSA{c}.roiArray{1}.metaData.subvolId)) ; continue ; end

            w1_t_pre = [w1_t_pre s_pre.caTSAArray.caTSA{c}.cellFeatures.get(w1_field)];
            w2_t_pre = [w2_t_pre s_pre.caTSAArray.caTSA{c}.cellFeatures.get(w2_field)];
        end
    else
        if (isempty( s_pre.caTSA.cellFeatures.get(w1_field))) ; s_pre.computeEvokedFeatures; end

        w1_t_pre = s_pre.caTSA.cellFeatures.get(w1_field);
        w2_t_pre = s_pre.caTSA.cellFeatures.get(w2_field);
    end 
    if (isfield(s_post.caTSAArray,'caTSA')) % are we dealing with mutliple subvolumes?
        w1_t_post = [];
        w2_t_post = [];
        for c=1:length(s_post.caTSAArray.caTSA) % this assumes pre/post are same size
            if (c == 1 && isempty( s_post.caTSAArray.caTSA{c}.cellFeatures.get(w1_field))) ; s_post.computeEvokedFeatures; end

            % matched subvolume?
            if (~isempty(subvol_id) & ~strcmp(subvol_id, s_post.caTSAArray.caTSA{c}.roiArray{1}.metaData.subvolId)) ; continue ; end

            w1_t_post = [w1_t_post s_post.caTSAArray.caTSA{c}.cellFeatures.get(w1_field)];
            w2_t_post = [w2_t_post s_post.caTSAArray.caTSA{c}.cellFeatures.get(w2_field)];            
        end
    else
        if (isempty( s_post.caTSA.cellFeatures.get(w1_field))) ; s_post.computeEvokedFeatures; end    

        w1_t_post = s_post.caTSA.cellFeatures.get(w1_field);
        w2_t_post = s_post.caTSA.cellFeatures.get(w2_field);       
    end

    % specific rejections
    invali = [];
    if (strcmp(ani_str, 'an014359'))
       invali = find(s_post.caTSAArray.caTSA{2}.ids >= 910040000); % an014359
       disp('an014359 ; dropping 091004');
    elseif (strcmp(ani_str, 'an017518'))
       invali = find(s_post.caTSAArray.caTSA{2}.ids >= 910040000); % an014359
       disp('an017518 ; dropping 091004');
    end

    % cells that meet criteria on either day are considered
    w1_pre_i = find(w1_t_pre >= score_thresh);
    w2_pre_i = find(w2_t_pre >= score_thresh);
    w1_post_i = find(w1_t_post >= score_thresh);
    w2_post_i = find(w2_t_post >= score_thresh);

    w1i = union(w1_pre_i,w1_post_i);
    w2i = union(w2_pre_i,w2_post_i);

    w1i = setdiff(w1i,invali);
    w2i = setdiff(w2i,invali);

    % show result ...
    mu_w1_pre = nanmean(w1_t_pre(w1i));
    mu_w2_pre = nanmean(w2_t_pre(w1i));
    mu_w1_post = nanmean(w1_t_post(w1i));
    mu_w2_post = nanmean(w2_t_post(w1i));

    disp(sprintf('w1 pre score: %0.3f post: %0.3f delta: %0.3f; w2 pre: %0.3f post: %0.3f delta: %0.3f', mu_w1_pre, mu_w1_post, (mu_w1_post-mu_w1_pre)/mu_w1_pre, ...
                                                                                                         mu_w2_pre, mu_w2_post, (mu_w2_post-mu_w2_pre)/mu_w2_pre ));

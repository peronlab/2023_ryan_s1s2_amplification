function plot_s1s2_projected_map (roi_array_wc, anim_data, dayi)
% plot_s1s2_projected_map('/Volumes/Data/s1s2_lesion/an016650/rois/an016650_2022_05_31_fov_001*.mat', all_dat.anim_data{find(strcmp(all_dat.anims, 'an016650'))}, 6)
    all_dat = get_s1s2_all_dat;
    fh = figure('Position', [0 0 400 400]);
    ax = axes('Position', [.05 .05 .9 .9]);
    
    % pull x/y/z from relevant roi array
    flist = dir(roi_array_wc);

    ids = [];
    roi_x = [];
    roi_y = [];
    
    roi_dir = fileparts_crossplatform(roi_array_wc);

    for f=1:length(flist)
        rA = load([roi_dir filesep flist(f).name]);
        rA = rA.obj;

        com = rA.getCOMs;
    
        roi_x = [roi_x com(:,1)'];
        roi_y = [roi_y com(:,2)'];
        ids = [ids rA.roiIds];
    end

    vali = find(ismember(anim_data.ids, ids));

    % now plot the maps
    w1_score = anim_data.day(dayi).w1ExclusiveTouchTrials;
    w2_score = anim_data.day(dayi).w2ExclusiveTouchTrials;
    dot_outside_range = 0;
    range_shown = [.05 .5];

    plot_projected_map (ax, roi_x, roi_y, w1_score(vali), w2_score(vali),  range_shown,  dot_outside_range, 0);
    ax.Title.String = [ax.Title.String ' ' strrep(anim_data.date_str{dayi},'_','-')];

function plot_projected_map (h, x, y, score_1, score_2, range_shown, dot_outside_range, show_border, ids)
% This will plot encoding map using circles of varying size and/or hue
%   h: handle to axes
%   x, y: vectors with coordinates for each cell ; z is ignored
%   score: vector with relevant parameter score for each cell 
%   range_shown: value range of score shown ; cells outside range are not shown; 2 el vector
%   dot_outside_range: if 1, cells outside range will be represented as dots
%   show_border: if 1, border shown ; default y; 2: border *only*
    if (nargin < 7) ; dot_outside_range = 1; end
    if (nargin < 8) ; show_border = 1; end
    if (nargin < 9) ; ids = [] ;end

    show_legend = 1;
    size_range = [5 500];

    % flip y
    y = -1*y;

    oscore_1 = score_1;
    oscore_1 = oscore_1-range_shown(1);
    oscore_1(find(oscore_1 < 0)) = 0;
    oscore_1 = oscore_1/diff(range_shown);
    oscore_1(find(oscore_1 > 1)) = 1;

    oscore_2 = score_2;
    oscore_2 = oscore_2-range_shown(1);
    oscore_2(find(oscore_2 < 0)) = 0;
    oscore_2 = oscore_2/diff(range_shown);
    oscore_2(find(oscore_2 > 1)) = 1;

    % eliminate baddies
    invali_1 = find(score_1 < range_shown(1) | score_1 > range_shown(2));
    score_1(invali_1) = nan;
    vali_1 = find(~isnan(score_1));

    invali_2 = find(score_2 < range_shown(1) | score_2 > range_shown(2));
    score_2(invali_2) = nan;
    vali_2 = find(~isnan(score_2));

    vali = union(vali_1,vali_2);
    invali = intersect(invali_1,invali_2);

    % rescale to [0 1] then size_range
    score_1 = score_1-range_shown(1);
    if (~isinf(range_shown(2)))
        score_1(find(score_1 > range_shown(2))) = range_shown(2);
    end
    score_1 = score_1/(max(score_1)/2);

    score_1 = score_1*(diff(size_range));
    score_1 = score_1+size_range(1);

    score_2 = score_2-range_shown(1);
    if (~isinf(range_shown(2)))
        score_2(find(score_2 > range_shown(2))) = range_shown(2);
    end
    score_2 = score_2/(max(score_2)/2);

    score_2 = score_2*(diff(size_range));
    score_2 = score_2+size_range(1);

    % and the plot ...
    hold(h,'on');
    
    % dotZ
    if (dot_outside_range)
        plot(h, x(invali), y(invali), 'k.', 'MarkerSize', 0.5);
    end

    % main action
    if (show_border ~= 2)
        cmat = [0.1+(oscore_1*.9) zeros(length(score_1),1) 0.1+(oscore_2*.9)];
        s = scatter(h, x(vali), y(vali), size_range(1)+diff(size_range)*((oscore_2(vali)+oscore_1(vali))/2), cmat(vali,:), 'fill');
    end
    if (show_border)
        s.MarkerFaceAlpha = 0.2;
        scatter(h, x(vali), y(vali), score(vali), color, 'LineWidth', 1);
    else
        s.MarkerFaceAlpha = 0.9;
    end
    if (show_border == 3)
        scatter(h, x(vali), y(vali), score(vali), color, 'LineWidth', 3);
    end

    % legend?
    if (show_legend)
       scatter(h, [20 20 20], -1*[50 100 150], [5 100 200], [0 0 0], 'fill');
       title(h, num2str(range_shown));
    end

    % ids?
    if (~isempty(ids))
        for c=1:length(ids)
            text(h, x(c), y(c), num2str(int64(ids(c))));
        end
    end

    axis(h, [0 512 -512 0]);
    aa = axis(h);
    set(h, 'TickDir', 'out', 'XTick', [aa(1) aa(2)], 'YTick', [aa(3) aa(4)]);    



function plot_2d_two_whiskers(all_dat, ani_str, di, ax, cell_type, touch, area_str, show_rejects, sat_val, size_range)
    if (nargin < 4 || isempty(ax)) ; figure; ax = axes ; end
    if (nargin < 5 || isempty(cell_type)) ; cell_type = 'usw' ;end
    if (nargin < 6 || isempty(touch)) ; touch = 'w1p' ;end
    if (nargin < 7 || isempty(area_str)) ; area_str = 's1_all' ;end
    if (nargin < 8 || isempty(show_rejects)) ; show_rejects = 0 ; end
    if (nargin < 9 || isempty(sat_val)) ; sat_val = 0.75 ; end
    if (nargin < 10 || isempty(size_range)) ; size_range = [1 25]; end
    
    % settings
    labels_shown = 1 ; 

    % prepare axes
    axes(ax);
    hold(ax,'on');

    % data
    ai = find(strcmp(ani_str, all_dat.anims));
    area_vali = get_s1s2_neuron_subset_idx(all_dat.anim_data{ai}.ids, area_str, all_dat, ai); % cell indices
    disp(sprintf('%d cells in %s tag', length(area_vali), area_str));

    vali_w1p = all_dat.types_by_idx{ai}.ev_touch_w1p_by_day{di};
    vali_w1r = all_dat.types_by_idx{ai}.ev_touch_w1r_by_day{di};
    vali_w2p = all_dat.types_by_idx{ai}.ev_touch_w2p_by_day{di};
    vali_w2r = all_dat.types_by_idx{ai}.ev_touch_w2r_by_day{di};

    vali_w1 = union(all_dat.types_by_idx{ai}.ev_touch_w1p_by_day{di}, all_dat.types_by_idx{ai}.ev_touch_w1r_by_day{di});
    vali_w2 = union(all_dat.types_by_idx{ai}.ev_touch_w2p_by_day{di}, all_dat.types_by_idx{ai}.ev_touch_w2r_by_day{di});

    score_valid_only_mat = nan*zeros(4,length(all_dat.anim_data{ai}.ids));
    if (1) % use mean dFF
        vali = all_dat.types_by_idx{ai}.touch_w1p_by_day{di};
        score_valid_only_mat(1,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_meanResponseTouchW1pExclusive(vali);
        vali = all_dat.types_by_idx{ai}.touch_w1r_by_day{di};
        score_valid_only_mat(2,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_meanResponseTouchW1rExclusive(vali);
        vali = all_dat.types_by_idx{ai}.touch_w2p_by_day{di};
        score_valid_only_mat(3,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_meanResponseTouchW2pExclusive(vali);
        vali = all_dat.types_by_idx{ai}.touch_w2r_by_day{di};
        score_valid_only_mat(4,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_meanResponseTouchW2rExclusive(vali);
        disp('mean')
    else % use P response
        vali = all_dat.types_by_idx{ai}.touch_w1p_by_day{di};
        score_valid_only_mat(1,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_probResponseTouchW1pExclusive(vali);
        vali = all_dat.types_by_idx{ai}.touch_w1r_by_day{di};
        score_valid_only_mat(2,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_probResponseTouchW1rExclusive(vali);
        vali = all_dat.types_by_idx{ai}.touch_w2p_by_day{di};
        score_valid_only_mat(3,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_probResponseTouchW2pExclusive(vali);
        vali = all_dat.types_by_idx{ai}.touch_w2r_by_day{di};
        score_valid_only_mat(4,vali) = all_dat.anim_data{ai}.day(di).evokedDffScores_probResponseTouchW2rExclusive(vali);
        disp('prob')
    end
    
    switch cell_type
        case 'usw'
            sat_color = all_dat.usw_color;
            ct_vali = all_dat.types_by_idx{ai}.usw_by_day{di};
        case 'bsw'
            sat_color = all_dat.bsw_color;
            ct_vali = all_dat.types_by_idx{ai}.bsw_by_day{di};

        case 'mw'
            sat_color = all_dat.mw_color;
            ct_vali = all_dat.types_by_idx{ai}.mw_by_day{di};

%         case 'usw'
%             sat_color = all_dat.usw_color;
%             ct_vali = all_dat.types_by_idx{ai}.usw_by_day{di};
    end

    switch touch
        case 'w1p'
            w_score = score_valid_only_mat(1,:);
        case 'w1r'
            w_score = score_valid_only_mat(2,:);
        case 'w2p'
            w_score = score_valid_only_mat(3,:);
        case 'w2r'
            w_score = score_valid_only_mat(4,:);
    end

    vali = intersect(area_vali, ct_vali);
    
    %w_score_debug= w_score;
    w_score = w_score/sat_val;
    w_score(find(w_score > 1)) = 1;
    w_score(find(w_score < 0)) = 0;
    
    xy = nan*zeros(length(w_score),2);
    xy(:,1) = all_dat.anim_data{ai}.x_um'; 
    xy(:,2) = all_dat.anim_data{ai}.y_um'; 
    %xyz(:,3) = all_dat.anim_data{ai}.z_um'; 

    % non-categorical
    if (show_rejects)
         invali = setdiff(1:length(xy), vali);
         invali = intersect(invali, area_vali);
         plot3(xyz(invali,1),xy(invali,2),xy(invali,3),'.','Color',  [1 1 1]*0.75, 'MarkerSize', 0.2);
     end

    disp(['Max depth: ' num2str(max(xy(:,3)))]);

    cmap_min = 0.75;
    cmap = [linspace(cmap_min*sat_color(1),sat_color(1),100); linspace(cmap_min*sat_color(2),sat_color(2),100); linspace(cmap_min*sat_color(3),sat_color(3),100)]';

    % categorical
    for v=1:length(vali)
        i = vali(v);

        w_color = cmap(max(1,round(w_score(i)*100)),:);
          
        col = w_color;
        col = min(col,[1 1 1]);
        
        msize(i) = round(max(size_range(1), (sqrt(w_score(i)))*size_range(2)));
        
         if w_score(i)>0
             plot3(xyz(i,1),xyz(i,2),xyz(i,3),'o','Color',col,'MarkerFaceColor',col,'MarkerSize',msize(i));
         end 
%         
%         if (show_rejects)
%         invali = setdiff(1:length(xyz), vali);
%         invali = intersect(invali, area_vali);
%         plot3(xyz(invali,1),xyz(invali,2),xyz(invali,3),'.','Color',  [1 1 1]*0.75, 'MarkerSize', 2);
%         end 
    end

    
    
    ax_w = 750;
    ax_h = 500%300; %500; %change depending on if 2 or 1 subvolumes imaged

    % pedestal
    pedestal_color = [1 1 1]*0.85;
    fill3([0 ax_w ax_w 0 ], [0 0 ax_w ax_w], [ 1 1 1 1 ]*500, pedestal_color, 'EdgeColor', pedestal_color);

    % touchups
    set(gca,'XTick',[],'YTick',[],'YDir','reverse','FontSize',15);
    %set(gca,'CameraPosition', [-3000 4000 -750]); % good for 1/2 as tall
    set(gca,'CameraPosition', [-3000 5000 0]);
    ax = gca;

    axis([0 ax_w 0 ax_w]);
    set(ax.YAxis,'visible','off');
    set(ax.XAxis,'visible','off');
    grid on  

    if (labels_shown)
        %ax.ZTick= [100 200 300 400];
        %ax.ZLabel.String = 'Depth \mum';
        ax.Title.String = sprintf('n=%d cells, %s %s', length(area_vali), ani_str, strrep(area_str,'_','-'));
    else
        %set(ax.ZAxis,'visible','off');
        %ax.ZTick= [];
    end
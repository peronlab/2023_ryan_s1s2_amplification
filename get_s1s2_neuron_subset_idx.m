function idx = get_s1s2_neuron_subset_idx(ids, restrict_tag, all_dat, ai)
% 
% Given a vector of neuron ids, it will return indices of ids that comply with the 
%  restrict tag.
%
%  restrict_tag: 's1_all','s1_manual_restrict', 's1_lesion'; same for s2
%
    
    switch restrict_tag
        case 's1_all'
            idx = find(ids >= all_dat.s1_id_range{ai}(1) & ids <= all_dat.s1_id_range{ai}(2));

        case 's1_lesion'
            idx = find(ids >= all_dat.s1_lesion_id_range{ai}(1) & ids <= all_dat.s1_lesion_id_range{ai}(2));
            
        case 'rad_effect'
            fname= [all_dat.root_dir filesep all_dat.anims{ai} filesep 'rois' filesep all_dat.anims{ai} '_lesion_roiIdList.mat' ];
            if (exist(fname, 'file'))
                idl = load(fname);
                idx = find(ismember(ids, idl.selectedRoiIds));
                idx = intersect(idx,find(ids >= all_dat.s1_id_range{ai}(1) & ids <= all_dat.s1_id_range{ai}(2)));
            else
                disp(['Could not find ' fname ' ; using range instead of file based restriction']);
                idx = find(ids >= all_dat.s1_id_range{ai}(1) & ids <= all_dat.s1_id_range{ai}(2));
            end
            
        case 's1_manual_restrict'
            fname = [all_dat.root_dir filesep all_dat.anims{ai} filesep 'rois' filesep all_dat.anims{ai} '_s1_roiIdList.mat' ];
            if (exist(fname, 'file'))
                idl = load(fname);
                idx = find(ismember(ids, idl.selectedRoiIds));
                idx = intersect(idx,find(ids >= all_dat.s1_id_range{ai}(1) & ids <= all_dat.s1_id_range{ai}(2)));
            else
                disp(['Could not find ' fname ' ; using range instead of file based restriction']);
                idx = find(ids >= all_dat.s1_id_range{ai}(1) & ids <= all_dat.s1_id_range{ai}(2));
            end
            %    idx = find(ids >= all_dat.s1_id_range{ai}(1) & ids <= all_dat.s1_id_range{ai}(2));

        case 's2_all'
            idx = find(ids >= all_dat.s2_id_range{ai}(1) & ids <= all_dat.s2_id_range{ai}(2));

        case 's2_lesion'
            idx = find(ids >= all_dat.s2_lesion_id_range{ai}(1) & ids <= all_dat.s2_lesion_id_range{ai}(2));
            
        case 's2_manual_restrict'
            fname = [all_dat.root_dir filesep all_dat.anims{ai} filesep 'rois' filesep all_dat.anims{ai} '_s2_roiIdList.mat' ];
            if (exist(fname, 'file'))
                idl = load(fname);
                idx = find(ismember(ids, idl.selectedRoiIds));
                idx = intersect(idx,find(ids >= all_dat.s2_id_range{ai}(1) & ids <= all_dat.s2_id_range{ai}(2)));
            else
                disp(['Could not find ' fname ' ; using range instead of file based restriction']);
                idx = find(ids >= all_dat.s2_id_range{ai}(1) & ids <= all_dat.s2_id_range{ai}(2));
            end
             %   idx = find(ids >= all_dat.s2_id_range{ai}(1) & ids <= all_dat.s2_id_range{ai}(2));
    end

        

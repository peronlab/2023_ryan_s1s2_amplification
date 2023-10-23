% script for evaluating redness threshold assignments - all it will do is color rois with redness > thresh purpl, all others blue
%
% dat = get_s1s2_single_animal_data;
% redness_eval(s, dat.ids, dat.day(3).red_score, 0.75)
%
%   Assumes session object "s" is loaded; bring up roi view and test this way. Use dat.day(x) and specify this in mdat.redness_day_idx.
%
function redness_eval (s, ids, redness, thresh)
    nr = 0;
    nn = 0;
    if (isfield(s.caTSAArray, 'caTSA'))
        for c=1:length(s.caTSAArray.caTSA)
            nr = nr+processSingleCaTSA(s.caTSAArray.caTSA{c}, ids, redness, thresh);
            nn = nn + s.caTSAArray.caTSA{c}.length;
        end
    else
        nr = processSingleCaTSA(s.caTSA, ids, redness, thresh);
        nn = nn + s.caTSA.length;
    end

    disp(sprintf('%0.3f of neurons are labeled', nr/nn));

function nr = processSingleCaTSA(caTSA, ids, redness, thresh)
    nr = 0;
    for rai=1:length(caTSA.roiArray)
        rA = caTSA.roiArray{rai};
        for r=1:length(rA.rois)
            rA.rois{r}.color = [0 0.5 1];
        end

        ii = find(ismember(ids, rA.roiIds));
        for i=1:length(ii)
            if (redness(ii(i)) > thresh)
                ri = find(rA.roiIds == ids(ii(i)));
                rA.rois{ri}.color = [1 0 1];
                nr = nr+1;
            end
        end
    end
    

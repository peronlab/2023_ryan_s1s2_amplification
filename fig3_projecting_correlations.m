% LR: march 22,2023
%wrapper for finding a subvolume with projection neurons and using it for
% correlations 
%functions used: find_usable_subvol (finds subvolume,aniamsl to discard. this is specific to this project), stim_locked_corrs_better_sub (correlations locked to stim),
%spont_corrs_better_sub(finds spont correlations), and
%corrs_allcellsametouch which gives a mean correlation of a MW cell to all
%other cells that respond to that same touch type! (projecting) 

%% Plot correlations locked to stimulus only in the subvolume with more cells for S1 img (S2 retro), projection cells only
[subs, cells, cells_subvol, discard]= find_usable_subvol(2,1); %retro 2 is s2 injection, proj only 1 means only use projection cells
S1_img_proj_bs= stim_locked_corrs_better_sub(2,1,subs,discard,0.5);
S1_img_nonproj_bs= stim_locked_corrs_better_sub(2,0,subs,discard, 0.5); % try non proj cells but use same subvols and such
%S1_img_proj_spont= spont_corrs_better_sub(2,1,subs,discard,0.8); 
%S1_img_proj_stimlock_allcell= corrs_allcellsametouch(2,1,subs,discard);
%% Plot correlations locked to stimulus only in the subvolume with more cells for S2 img (S1 retro), projection cells only
% not removing bidir cells here but won't include in anything!! 
[subs1, cells1, cells_subvol1, discard1]= find_usable_subvol(1,1); 
S2_img_proj_bs= stim_locked_corrs_better_sub(1,1,subs1,discard1,0.5);
S2_img_nonproj_bs= stim_locked_corrs_better_sub(1,0,subs1,discard1,0.5); % try non proj cells but use same subvols and such
%S2_img_proj_spont= spont_corrs_better_sub(1,1,subs1,discard1,0.8); 
%S2_img_proj_stimlock_allcell= corrs_allcellsametouch(1,1,subs1,discard1);
%% nonprojecting
%S1_img_non_proj_bs= stim_locked_corrs_better_sub(2,0,subs,discard, );
%S1_img_non_proj_spont= spont_corrs_better_sub(2,0,subs,discard,0.8); 
%S2_img_non_proj_bs= stim_locked_corrs_better_sub(1,0,subs1,discard1);
%S2_img_non_proj_spont= spont_corrs_better_sub(1,0,subs1,discard1,0.8); 
%% stats testing
%% 
%S1_corrs (:,:,1)= S1_img_proj_bs; %stim locked correlations w/in groups, projection only neurons
S1_corrs (:,:,1)= S1_img_proj_spont;%spontaenous correlaions w/in groups, projection only
%S1_corrs (:,:,3)= S1_img_proj_stimlock_allcell;%correlations of projection cells- for each group, it is mean correlation of that cell type with all cells that respond to the same touch type

S2_corrs (:,:,1)= S2_img_proj_bs; %stim locked correlations w/in groups, projection only neurons
S2_corrs (:,:,2)= S2_img_proj_spont;%spontaenous correlaions w/in groups, projection only
%S2_corrs (:,:,3)= S2_img_proj_stimlock_allcell;%correlations of projection cells- for each group, it is mean correlation of that cell type with all cells that respond to the same touch type

for i=1:2 
    if i==1
        x= S1_corrs;
    else 
        x= S2_corrs; 
    end 
    for j=1:2%3 %loop through the types of correlatrions 
        A= x(1,:,j);
        B= x(2,:,j);
        C= x(3,:,j);
       [h, pvab]= ttest(A,B);
       [h, pvac]= ttest(A,C);
       [h, pvbc]= ttest(B,C);
       pvals (j,:,i) = [pvab, pvac, pvbc];%ab, ac, bc
    end 
end 
% A=S2_img_proj_stimlock_allcell(1,:);
% B=S2_img_proj_stimlock_allcell(2,:);
% C=S2_img_proj_stimlock_allcell(3,:);
% [h, pvab]= ttest(A,B);
% [h, pvac]= ttest(A,C);
% [h, pvab]= ttest(B,C);
%% 
% run this script to automatically adjust whisker gui params for this
% project -- 2whisker coarse resolution. 
 %% 1. first you will need to go into the animal whisker folder, and run this from the wta folder. 
info=dir ('*wta.mat') % lists all the wta.mat files. we will run the adjust for each of them 

%%
for i= 1:length (info) %for each date 
    load (info(i).name) % what to load 

    o.params.detectContactParams.kappaThreshMult = [5 1]; %adjustments for contact detection for this particular project 
    o.params.detectContactParams.doWhisker = [0 1];
    o.detectContacts;
% if doing general adjustments, comment out these 3 commands-- dont need,
% only for this 2 whisker project. 

    o.classifyContactsProRet; %then run this as a quick fix-- usually helps (in general case)
    o.manuallySetBarFracInReach(0.8); % another good adjustment in general 

    o.saveToFile; 
    %clear all; 
end 

%% this isn't going to be automated. run individaully to open each gui 
% for each, reload and then: 
% o.guiTrialBrowser;
% o.guiData.trialBrowser.maxDToBar = 20;
% o.updatePaths('~/Desktop/tmp_whisker/an014359/2022_08_25_vidwhisk');
% at the end, save, reupload wta folder to cluster.
%% Background

% Naureen Ghani gave story about mice doing wiggles for low-contrast
% stimuli, she's trying to publish just a behavioral story but would like
% neural data if possible. I was checking if something might be here in the
% widefield for her to use. Slightly different protocol (she looked at IBL
% data which doesn't have interactive delay), but if anything that might
% make this dataset cleaner. e.g. was looking at whether incidence of
% wiggles at interactive onset only happens on trials with lower initial
% visual cortex activity.


%% Load data
load("C:\Users\petersa\Documents\Previous_labs\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\paper_unused\trial_activity_choiceworld_wfonly.mat");
load("C:\Users\petersa\Documents\Previous_labs\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat");

%% Look at wheel wiggles and widefield

% Grab data
t = trial_data_all.t;

x = cat(1,trial_data_all.trial_info_all{:});
stim_cat = cell2mat(cellfun(@(x) x.stimulus,x,'uni',false));

x = cat(1,trial_data_all.wheel_all{:});
wheel_cat = cell2mat(x);

n_vs = 200;
wf_cat = AP_deconv_wf(cell2mat(cat(1,trial_data_all.fluor_all{:})));
wf_cat = wf_cat - nanmean(wf_cat(:,t<0,:),2);


% define wiggles as wheel velocity reversals
rev = sign(wheel_cat(:,1:end-1)) ~= sign(wheel_cat(:,2:end));
rev_grp = ap.groupfun(@mean,rev,stim_cat,[]);

figure;

subplot(1,2,1);
plot(t(2:end),rev_grp')
colororder(gca,ap.colormap('BKR',11));
xlabel('Time from stim');
ylabel('Frac trials with reversals')
xline([0,0.2],'g');
xline([0.5,0.7],'m');

subplot(1,2,2); hold on;
plot(unique(stim_cat),nanmean(rev_grp(:,t > 0.5 & t < 0.7),2),'m');
plot(unique(stim_cat),nanmean(rev_grp(:,t > 0 & t < 0.2),2),'g');
xlabel('Contrast');
ylabel('Reversals');

% (flip wf to combine neg trials)
% (used to use AP_reclect_widefield - align file is elsewhere now)
wf_cat_mirror = wf_cat;

right_stim_trials = stim_cat < 0;
mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(fliplr(U_master(:,:,1:n_vs)),[],n_vs);
wf_cat_mirror(right_stim_trials,:,:) = reshape(transpose( ...
    mirror_matrix*reshape(wf_cat(right_stim_trials,:,:),[],n_vs)'), ...
    size(wf_cat(right_stim_trials,:,:)));


% Split by reversals or not
use_trials = find(abs(stim_cat) == 0);
x_grp = any(rev(use_trials,t>0.5 & t<0.7),2);

a = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(ap.groupfun(@nanmean,wf_cat_mirror(use_trials,:,:),x_grp,[],[]),[3,2,1]));
ap.imscroll(a,t);axis image;
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));

wheel_grp = ap.groupfun(@nanmean,wheel_cat(use_trials,:),x_grp,[]);
rev_grp = ap.groupfun(@nanmean,rev(use_trials,:),x_grp,[]);
figure;
subplot(1,2,1);plot(t,wheel_grp');
subplot(1,2,2);plot(t(2:end),rev_grp');



% split trials by reversal amount
use_trials = find(stim_cat == -0.06);
x = nanmean(rev(use_trials,t>0.5 & t<0.7),2);

n_split = 2;
x_grp = discretize(tiedrank(x)/length(x),linspace(0,1,n_split+1));

a = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(ap.groupfun(@mean,wf_cat(use_trials,:,:),x_grp,[],[]),[3,2,1]));
ap.imscroll(a,t);axis image;
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));

wheel_grp = ap.groupfun(@nanmean,wheel_cat(use_trials,:),x_grp,[]);
rev_grp = ap.groupfun(@nanmean,rev(use_trials,:),x_grp,[]);
figure;
subplot(1,2,1);plot(t,wheel_grp');
subplot(1,2,2);plot(t(2:end),rev_grp');


% potentially looks like: little wiggles = big initial response and small
% response on beep, more wiggles = small initial response and bigger
% response on beep


% split by vis activity and look at wiggle incidence

% (plot avg wf by stim)
a = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(ap.groupfun(@nanmean,wf_cat,stim_cat,[],[]),[3,2,1]));
ap.imscroll(a,t);axis image;
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));

% (plot avg wf by stim - combine sides)
a = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(ap.groupfun(@nanmean,wf_cat_mirror,abs(stim_cat),[],[]),[3,2,1]));
ap.imscroll(a,t);axis image;
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));



% (after drawing ROI, grab fluorescence)
roi_trace = permute(ap.wf_roi(U_master(:,:,1:n_vs), ...
    permute(wf_cat,[3,2,1]),[],[],roi.mask),[3,2,1]);

% (split data within animal)
n_split = 4;
use_t = t >= 0.05 & t <= 0.15;
use_trials = find(stim_cat == 0.06);

trial_grp = discretize(tiedrank(nanmean(roi_trace(use_trials,use_t),2))/length(use_trials),linspace(0,1,n_split+1));

a = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(ap.groupfun(@nanmean,wf_cat(use_trials,:,:),trial_grp,[],[]),[3,2,1]));
ap.imscroll(a,t);axis image;
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));

rev_grp = ap.groupfun(@nanmean,rev(use_trials,:),trial_grp,[]);
figure;plot(t(3:end),rev_grp')

wheel_grp = ap.groupfun(@nanmean,wheel_cat(use_trials,:),trial_grp,[]);










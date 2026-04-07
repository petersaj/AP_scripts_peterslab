% Cleaning up PG pupil code to be consistent with other figure code


%% Load general data

%%% Load data for figure
load_dataset = 'noact';
Marica_2026.figures.load_data;
%%%

%% Load pupil data
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2026','data');
pupil_data_filename = fullfile(data_path,'pupil_passive');
load(pupil_data_filename);

mousecam_framerate = 30;
pupil_t = -0.5:1/mousecam_framerate:1.5;
pupil_deriv_t = pupil_t(1:end-1) + diff(pupil_t)/2;

plot_day_bins = [-Inf,-2,0,2,Inf];
pupil_plot_day_grp = discretize(max(pupil_diameter_grp.ld,-inf),plot_day_bins);

% Concatenate pupil data and create indicies
pupil_diameter = cell2mat(pupil.pupil_diameter);
% (fully nan-out any trials with missing data)
pupil_diameter = pupil_diameter.*ap.nanout(any(isnan(pupil_diameter),2));
% (get pupil diameter derivative)
pupil_diameter_deriv = diff(pupil_diameter,[],2);

pupil_diameter_grp.animal = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(grp2idx(pupil.animal)),pupil.pupil_diameter,'uni',false));

pupil_diameter_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(bhv.days_from_learning),pupil.pupil_diameter,'uni',false));

pupil_diameter_grp.stim = cell2mat(pupil.trial_stim_values);

% Plot average pupil diameter derivative
[pupil_deriv_avg,pupil_deriv_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},pupil_diameter_deriv, ...
    pupil_diameter_grp.animal,[pupil_plot_day_grp,pupil_diameter_grp.stim]);

pupil_deriv_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},pupil_diameter_deriv, ...
    pupil_diameter_grp.animal,[pupil_plot_day_grp,pupil_diameter_grp.stim]);

stim_color = ap.colormap('BKR',3);
unique_stim = unique(pupil_diameter_grp.stim);

figure;
h = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
for curr_day = 1:length(plot_day_bins)-1
    nexttile; hold on;
    axis square;

    for curr_stim = unique_stim'
        curr_data_idx = ismember(pupil_deriv_avg_grp,[curr_day,curr_stim],'rows');
        ap.errorfill(pupil_deriv_t,pupil_deriv_avg(curr_data_idx,:), ...
            pupil_deriv_sem(curr_data_idx,:),stim_color(unique_stim==curr_stim,:), ...
            0.5,true,2);
    end

end

linkaxes(h.Children,'xy');
ylim([-0.03, 0.035]);
ap.prettyfig;

% Plot pupil derivative AUC
pupil_auc_t = [0,1];
use_pupil_t = isbetween(pupil_deriv_t,pupil_auc_t(1),pupil_auc_t(2));
pupil_deriv_AUC = trapz(pupil_deriv_t(use_pupil_t),pupil_diameter_deriv(:,use_pupil_t),2);

[pupil_deriv_AUC_avg,pupil_deriv_AUC_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},pupil_deriv_AUC, ...
    pupil_diameter_grp.animal,[pupil_plot_day_grp,pupil_diameter_grp.stim]);

pupil_deriv_AUC_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},pupil_deriv_AUC, ...
    pupil_diameter_grp.animal,[pupil_plot_day_grp,pupil_diameter_grp.stim]);

% Plot across day bins
figure; hold on
set(gca,'ColorOrder',stim_color)
for curr_stim = unique_stim'
    curr_data_idx = ismember(pupil_deriv_AUC_grp(:,2),curr_stim);
    errorbar([],pupil_deriv_AUC_avg(curr_data_idx), ...
        pupil_deriv_AUC_sem(curr_data_idx),'LineWidth',2,'CapSize',0)
end
axis square padded;
ylabel('Mean derivative AUC')
 

% ~~~ STATS ~~~
% (compare day i to i+1)
print_stat('\n--FIG S8--\n');
print_stat('Pupil deriv max\n');

[pupil_deriv_AUC_dayavg,pupil_deriv_AUC_dayavg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},pupil_deriv_AUC, ...
    (1:size(pupil_deriv_AUC,1))', ...
    [pupil_diameter_grp.animal,pupil_plot_day_grp,pupil_diameter_grp.stim]);

for curr_compare_day = 1:length(plot_day_bins)-2

    compare_day_grps = curr_compare_day+[0,1];

    pupil_stat_usedata = ismember(pupil_deriv_AUC_dayavg_grp(:,2),compare_day_grps);
    [pupil_stat_meas,pupil_stat_grp] = ap.nestgroupfun({@nanmean,@diff}, ...
        pupil_deriv_AUC_dayavg(pupil_stat_usedata,:,:), ...
        pupil_deriv_AUC_dayavg_grp(pupil_stat_usedata,2), ...
        pupil_deriv_AUC_dayavg_grp(pupil_stat_usedata,3));

    [~,~,pupil_shuff_grp] = unique(pupil_deriv_AUC_dayavg_grp(pupil_stat_usedata,[1,3]),'rows');

    n_shuff = 10000;
    pupil_stat_null = nan(size(pupil_stat_meas,1),n_shuff);
    for curr_shuff = 1:n_shuff
        pupil_data_shuff = ap.shake(pupil_deriv_AUC_dayavg(pupil_stat_usedata,:,:),1,pupil_shuff_grp);
        pupil_stat_null(:,curr_shuff) = ap.nestgroupfun({@nanmean,@diff},pupil_data_shuff, ...
            pupil_deriv_AUC_dayavg_grp(pupil_stat_usedata,2),pupil_deriv_AUC_dayavg_grp(pupil_stat_usedata,3));
    end

    pupil_stat_rank = permute(tiedrank(permute([pupil_stat_meas,pupil_stat_null],[2,1,3])),[2,1,3]);
    pupil_stat_p = 1-pupil_stat_rank(:,1,:)/(n_shuff+1);

    print_stat('pupil: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(pupil_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_stim = unique(pupil_stat_grp(:,1))'
        curr_stat_idx = ismember(pupil_stat_grp,curr_stim,'rows');
        print_stat('Pupil, Stim %3.f, p = %.2g%s\n', ...
            pupil_stat_grp(curr_stat_idx),pupil_stat_p(curr_stat_idx),stat_sig(curr_stat_idx));
    end
end

%%%%% WORKING HERE
% Stim 0 isn't sig pre-learning? that's a little surprising




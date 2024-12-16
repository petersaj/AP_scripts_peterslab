%% ~~ Grab and package trial data
%
% Structure 'trial_data' is filled in each cell

%% Set parameters

% Trial time to grab
trial_t_window = [-0.5,1];
trial_t_rate = 50;
trial_t = trial_t_window(1):1/trial_t_rate:trial_t_window(2);

% Widefield components to keep 
n_wf_components = 200;

% Striatum multiunit chunk size
striatum_mua_length = 200; % um

% Check if task or passive
if contains(bonsai_workflow,'wheel')
    task_workflow = true;
elseif contains(bonsai_workflow,'passive')
    task_workflow = false;
else
    error('Bonsai workflow "%s": not "wheel" or "passive"',bonsai_workflow);
end

% Create structure for trial data
trial_data = struct;

%% Trial information

if task_workflow
    % Task
    align_t = reshape(stimOn_times(1:n_trials),[],1);
    trial_data.stim_to_move = stim_to_move(1:n_trials);
    trial_data.outcome = vertcat(trial_events.values(1:n_trials).Outcome);

elseif ~task_workflow
    % Passive
    n_stim = length(stimOn_times);

    align_t = reshape(stimOn_times,[],1);
    stim_x_cat = vertcat(trial_events.values.TrialStimX);
    trial_data.stim = stim_x_cat(1:n_stim); 
end

% Set alignment sample times
interp_t = align_t + trial_t;

%% Wheel

trial_data.wheel_velocity = interp1(timelite.timestamps, ...
    wheel_velocity,interp_t);

trial_data.wheel_move = interp1(timelite.timestamps, ...
    +wheel_move,interp_t,'previous');


%% Striatum MUA (in length chunks)

% Find striatum boundaries
AP_longstriatum_find_striatum_depth

if ~any(isnan(striatum_depth))

% Discretize spikes by depth
depth_group_edges = striatum_start:striatum_mua_length:striatum_end;
spike_depth_group = discretize(spike_depths,depth_group_edges);

n_depths = max(spike_depth_group);

% Get trial PSTH
[~,trial_data.striatum_mua] =  ap.psth(spike_times_timelite, ...
    align_t,spike_depth_group, ...
    'window',trial_t_window,'bin_size',1/trial_t_rate);

end

%% Cortical widefield

trial_data.widefield = interp1(wf_t,wf_V(1:n_wf_components,:)',interp_t,'previous');




























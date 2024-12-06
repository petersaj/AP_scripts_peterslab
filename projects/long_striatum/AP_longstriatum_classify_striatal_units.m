% Classify striatal cell types
% Only for quality metric "single units": MSN, FSI, TAN

% (Should be done elsewhere - this finds striatal depth on probe)
if ~exist('striatum_depth','var')
    AP_longstriatum_find_striatum_depth
end

% Only use quality-controlled "single units", within striatum
striatal_single_units = ...
    strcmp(template_qc_labels,'singleunit') & ...
    ~isnan(discretize(template_depths,striatum_depth));

if verbose; fprintf('Ephys: Classifying striatal cells...\n'); end

% Get spike acgs (messy for now - hard-code initializing and running only
% for striatal, since this takes a little while)
spike_acg = nan(size(templates,1),2001);
spike_acg(striatal_single_units,:) = cell2mat(arrayfun(@(x) ...
    ap.ephys_spike_cg(x),find(striatal_single_units),'uni',false));

% Get time to get to 90% steady-state value
acg_steadystate = nan(size(templates,1),1);
acg_steadystate(striatal_single_units) = arrayfun(@(x) ...
    find(spike_acg(x,ceil(size(spike_acg,2)/2):end) > ...
    mean(spike_acg(x,end-100:end),2)*0.9,1,'first'),find(striatal_single_units));

% Get average firing rate from whole session
spike_rate = accumarray(spike_templates,1)/diff(prctile(spike_times_timelite,[0,100]));

% Define cell types
% (NOTE: Julie uses acg_steadystate of 40, seemed better here for 20)
striatum_celltypes = struct;

striatum_celltypes.msn = striatal_single_units & ... % striatal single unit
    waveform_duration_peaktrough >= 400 & ... wide waveform
    acg_steadystate < 20; % fast time to steady state

striatum_celltypes.fsi = striatal_single_units & ... % striatal single unit
    waveform_duration_peaktrough < 400 & ... % narrow waveform
    acg_steadystate < 20; % slow time to steady state

% striatum_celltypes.tan = striatal_single_units & ... % striatal single unit
%     spike_rate >= 4 & spike_rate <= 12 & ... % steady firing rate
%     waveform_duration_peaktrough >= 400 & ... wide waveform
%     acg_steadystate >= 20; % slow time to steady state

% !! NOT USING WAVEFORM DURATION HERE - some clear TANs with short wfs
striatum_celltypes.tan = striatal_single_units & ... % striatal single unit
    spike_rate >= 4 & spike_rate <= 12 & ... % steady firing rate
    acg_steadystate >= 20; % slow time to steady state

% Plot units by cell type
if false

    str_celltype_colors = lines(4);
    str_unit_colors = str_celltype_colors( ...
        sum([striatum_celltypes.msn, ...
        striatum_celltypes.fsi, ...
        striatum_celltypes.tan].*[1:3],2)+1,:);

    figure; tiledlayout(1,2);

    % Scatter plot units by properties
    nexttile; hold on
    scatter3( ...
        waveform_duration_peaktrough(striatal_single_units), ...
        acg_steadystate(striatal_single_units), ...
        spike_rate(striatal_single_units), ...
        20,str_unit_colors(striatal_single_units,:),'filled');
    xlabel('Peak-trough duration');
    ylabel('ACG steady state time');
    zlabel('Spike rate');
    set(gca,'xscale','log','zscale','log');
    view(45,45);
    axis tight vis3d;

    % Plot average ACG by cell type
    str_grp_idx = sum([striatum_celltypes.msn, ...
        striatum_celltypes.fsi, ...
        striatum_celltypes.tan].*[1:3],2);

    spike_acg_grp = ap.groupfun(@mean,normalize(spike_acg,2,'range'),str_grp_idx,[]);

    nexttile;plot(spike_acg_grp(2:end,:)');
    ylabel('Range-normalized ACG');
    legend({'MSN','FSI','TAN'});

end













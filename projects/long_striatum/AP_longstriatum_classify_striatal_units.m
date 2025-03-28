% Classify striatal cell types
% Only for quality metric "single units": MSN, FSI, TAN

% (Should be done elsewhere - this finds striatal depth on probe)
if ~exist('striatum_depth','var')
    AP_longstriatum_find_striatum_depth
end

% Only use quality-controlled "single units", within striatum
% (REMOVED FROM BELOW - this can be done at later step)
striatal_units = ~isnan(discretize(template_depths,striatum_depth));
single_units = ismember(template_qc_labels,{'singleunit'});

if verbose; fprintf('Ephys: Classifying striatal cells...\n'); end

% %% Attempt at metrics w/o bombcell
% 
% %%%% USE RAW WAVEFORMS?
% % 
% % raw_templates_fn = fullfile(qMetrics_path,'templates._bc_rawWaveforms.npy');
% % raw_templates = permute(readNPY(raw_templates_fn),[1,3,2]);
% % 
% % [~,raw_max_site] = max(max(abs(raw_templates),[],2),[],3);
% % raw_waveforms = cell2mat(arrayfun(@(x) raw_templates(x,:,raw_max_site(x)), ...
% %     (1:size(raw_templates,1))','uni',false));
% % 
% % raw_waveforms = raw_waveforms(good_templates,:);
% % 
% %%%

% Get spike acgs (messy for now - hard-code initializing and running only
% for striatal, since this takes a little while)
spike_acg = nan(size(templates,1),2001);
spike_acg(striatal_units,:) = cell2mat(arrayfun(@(x) ...
    ap.ephys_spike_cg(x),find(striatal_units),'uni',false));

% Get time to get to 90% steady-state value
acg_steadystate = nan(size(templates,1),1);
acg_steadystate(~any(isnan(spike_acg),2)) = arrayfun(@(x) ...
    find(spike_acg(x,ceil(size(spike_acg,2)/2):end) > ...
    mean(spike_acg(x,end-100:end),2)*0.9,1,'first'),find(~any(isnan(spike_acg),2)));

% (UNUSED: ACG RATIO)
acg_early = max(spike_acg(:,1001:1001+300),[],2);
acg_late = max(spike_acg(:,end-200:end-100),[],2);
acg_ratio = acg_late./acg_early;

% Get average firing rate from whole session
spike_rate = accumarray(spike_templates,1)/diff(prctile(spike_times_timelite,[0,100]));

% Define cell types
% (NOTE: Julie uses acg_steadystate of 40, seemed better here for 20)
striatum_celltypes = struct;

striatum_celltypes.msn = ... % striatal_single_units & ... % striatal single unit
    waveform_duration_peaktrough >= 400 & ... wide waveform
    acg_steadystate < 20; % fast time to steady state

striatum_celltypes.fsi = ... % striatal_single_units & ... % striatal single unit
    waveform_duration_peaktrough < 400 & ... % narrow waveform
    acg_steadystate < 20; % slow time to steady state

% striatum_celltypes.tan = striatal_single_units & ... % striatal single unit
%     spike_rate >= 4 & spike_rate <= 12 & ... % steady firing rate
%     waveform_duration_peaktrough >= 400 & ... wide waveform
%     acg_steadystate >= 20; % slow time to steady state

% !! NOT USING WAVEFORM DURATION HERE - some clear TANs with short wfs
striatum_celltypes.tan = ... % striatal_single_units & ... % striatal single unit
    spike_rate >= 4 & spike_rate <= 16 & ... % steady firing rate
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
        waveform_duration_peaktrough(striatal_units), ...
        acg_steadystate(striatal_units), ...
        spike_rate(striatal_units), ...
        20,str_unit_colors(striatal_units,:),'filled');
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


%% Julie's method
% EPHYS PROPERTIES BROKEN AT THE MOMENT - BUG?
% 
% ephysProperties = bc.ep.runAllEphysProperties(convertStringsToChars(kilosort_path), convertStringsToChars(qMetrics_path), false,[]);
% ephysProperties = ephysProperties(good_templates,:);
% 
% striatum_celltypes = struct;
% 
% spike_param_wide_narrow_wavelength = 400;
% spike_param_tan_min_wavelength = 300;
% spike_param_post_spike_suppression = 30;
% spike_param_prop_isi = 0.1;
% 
% striatum_celltypes.msn = striatal_single_units & ...
%     ephysProperties.waveformDuration_peakTrough_us > spike_param_wide_narrow_wavelength &...
%     ephysProperties.postSpikeSuppression_ms < spike_param_post_spike_suppression;
% 
% striatum_celltypes.fsi = striatal_single_units & ...
%     ephysProperties.waveformDuration_peakTrough_us <= spike_param_wide_narrow_wavelength &...
%     ephysProperties.propLongISI <= spike_param_prop_isi;
% 
% striatum_celltypes.tan = striatal_single_units & ...
%     ephysProperties.waveformDuration_peakTrough_us > spike_param_wide_narrow_wavelength &...
%     ephysProperties.postSpikeSuppression_ms >= spike_param_post_spike_suppression;
% 
% striatum_celltypes.uin = striatal_single_units & ...
%     ephysProperties.waveformDuration_peakTrough_us <= spike_param_wide_narrow_wavelength &...
%     ephysProperties.propLongISI > spike_param_prop_isi;
% 
% % % SLIGHTLY CHANGED: remove wide TAN requirement (some clear TANs with
% % % narrow ~300 wave duration?)
% % 
% % striatum_celltypes_bombcell.msn = striatal_single_units & ...
% %     ephysProperties.waveformDuration_peakTrough_us > spike_param_wide_narrow_wavelength &...
% %     ephysProperties.postSpikeSuppression_ms < spike_param_post_spike_suppression;
% % 
% % striatum_celltypes_bombcell.fsi = striatal_single_units & ...
% %     ephysProperties.waveformDuration_peakTrough_us <= spike_param_wide_narrow_wavelength &...
% %     ephysProperties.postSpikeSuppression_ms < spike_param_post_spike_suppression;
% % 
% % striatum_celltypes_bombcell.tan = striatal_single_units & ...
% %     ephysProperties.postSpikeSuppression_ms >= spike_param_post_spike_suppression & ...
% %     ephysProperties.waveformDuration_peakTrough_us >= spike_param_tan_min_wavelength & ...
% %     ~isnan(discretize(ephysProperties.mean_firingRate,[4,16]));
% 
% 
% % Plot units by cell type
% if false
% 
%     str_celltype_colors = lines(4);
%     str_unit_colors = str_celltype_colors( ...
%         sum([striatum_celltypes.msn, ...
%         striatum_celltypes.fsi, ...
%         striatum_celltypes.tan].*[1:3],2)+1,:);
% 
%     % Scatter plot units by properties
%     figure; hold on
%     scatter3( ...
%         ephysProperties.waveformDuration_peakTrough_us(striatal_single_units), ...
%         ephysProperties.postSpikeSuppression_ms(striatal_single_units), ...
%         ephysProperties.mean_firingRate(striatal_single_units), ...
%         20,str_unit_colors(striatal_single_units,:),'filled');
%     xlabel('Peak-trough duration');
%     ylabel('Post-spike suppression');
%     zlabel('Firing rate');
% 
%     view(45,45);
%     axis tight vis3d
% 
% end











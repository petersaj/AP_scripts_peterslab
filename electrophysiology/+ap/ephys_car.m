function ephys_car(data_raw_filenames,data_car_filename,probe_info)
% ephys_car(data_raw_filename,data_car_filename,probe_info)
% Apply common average referencing to Neuroipxels data
%
% 1) Get median across simultaneously-sampled channels (1/bank)
% 2) Get scaling factor of median for each channel
% 3) Scale and subtract median for each channel
%
% probe_info: output from plab.ephys.oe_probe_info


% Hardcode channel number and get channel index across ADCs
% (X ADCs, Y channels each, same-index channels sampled together e.g.
% ADC1/Ch1 + ADC2/Ch1...)
% (https://open-ephys.github.io/gui-docs/User-Manual/Plugins/Neuropixels-CAR.html)

n_chan = probe_info.n_chan;

probe_version = sscanf(probe_info.probe_type,'%*s%d%*s');
switch probe_version
    case 1
        % Neuropixels 1.0
        n_adcs = 32;
    case 2
        % Neuropixels 2.0
        n_adcs = 24;
end

ch_adc_idx = reshape(reshape(1:n_chan,2,[])',n_chan/n_adcs,[]);

% Choose data chunk size
chunkSize = 1e6;

% Open CAR file for writing
fid_car = fopen(data_car_filename, 'w');

for curr_recording = 1:length(data_raw_filenames)

    curr_data_raw_filename = data_raw_filenames(curr_recording);
    fprintf('CAR-ing data: %s (file %d/%d) \n',curr_data_raw_filename, ...
        curr_recording,length(data_raw_filenames))

    % Get number of chunks in file (for printing progress)
    ap_filename_dir = dir(curr_data_raw_filename);
    nSampsTotal = ap_filename_dir.bytes/n_chan/2;
    nChunksTotal = ceil(nSampsTotal/chunkSize);

    % Open raw file for reading
    fid_raw = fopen(curr_data_raw_filename, 'r');

    % Loop through data chunks and apply common average referencing
    chunkInd = 1;
    data_eof = feof(fid_raw);
    while ~data_eof

        % Load in current data chunk
        curr_data = fread(fid_raw, [n_chan chunkSize], '*int16');

        % Subtract median across channels
        curr_data_centered = curr_data - median(curr_data,2);

        % Loop through ADC channel indicies and shanks
        curr_data_car = zeros(size(curr_data),'int16');
        for curr_adc_idx = 1:size(ch_adc_idx,1)

            %%%% IN PROGRESS: BY SHANK?
            % for curr_shank = unique(probe_info.kcoords)'
            %     % Get same-index channels across ADCs within shank
            %     curr_channels = intersect(ch_adc_idx(curr_adc_idx,:), ...
            %         find(probe_info.kcoords == curr_shank));

            curr_channels = ch_adc_idx(curr_adc_idx,:);

            % Get median signal across channels
            curr_data_adc_median = median(curr_data_centered(curr_channels,:),1);

            % Get scaling factor of median for each channel
            curr_data_adc_median_scale = ...
                sum(curr_data_centered(curr_channels,:).*curr_data_adc_median,2)./ ...
                sum(curr_data_adc_median.^2);
            curr_data_adc_median_scaled = int16(curr_data_adc_median_scale.* ...
                double(curr_data_adc_median));

            % Subtract scaled median from data
            curr_data_car(curr_channels,:) = curr_data_centered(curr_channels,:) - ...
                curr_data_adc_median_scaled;
            % end
        end

        fwrite(fid_car, curr_data_car, 'int16');

        ap.print_progress_fraction(chunkInd,nChunksTotal);

        chunkInd = chunkInd+1;
        data_eof = feof(fid_raw);

    end

    % Close raw file
    fclose(fid_raw);

end

% Close CAR file
fclose(fid_car);

disp('Finished CAR');





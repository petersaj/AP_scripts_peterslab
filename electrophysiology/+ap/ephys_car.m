function ephys_car(data_raw_filename,data_car_filename)
% ephys_car(data_raw_filename,data_car_filename)
% Apply common average referencing to Neuroipxels data
%
% 1) Get median across simultaneously-sampled channels (1/bank)
% 2) Get scaling factor of median for each channel
% 3) Scale and subtract median for each channel

% Allow for multiple filenames: package in cell if only one filename
if ~iscell(data_raw_filename)
    data_raw_filename = {data_raw_filename};
end

% Hardcode channel number and get channel index across ADCs
% (16 ADCs, 12 channels each, same-index channels sampled together e.g.
% ADC1/Ch1 + ADC2/Ch1...)
% (https://open-ephys.github.io/gui-docs/User-Manual/Plugins/Neuropixels-CAR.html)
n_chan = 384;
ch_adc_idx = reshape(reshape(1:n_chan,2,[])',12,[]);

% Choose data chunk size
chunkSize = 1e6;

% Open CAR file for writing
fid_car = fopen(data_car_filename, 'w');

for curr_recording = 1:length(data_raw_filename)

    curr_data_raw_filename = data_raw_filename{curr_recording};
    fprintf('CAR-ing data: %s (file %d/%d) \n',curr_data_raw_filename, ...
        curr_recording,length(data_raw_filename))

    % Get number of chunks in file (for printing progress)
    ap_filename_dir = dir(curr_data_raw_filename);
    nSampsTotal = ap_filename_dir.bytes/n_chan/2;
    nChunksTotal = ceil(nSampsTotal/chunkSize);

    % Open raw file for reading
    fid_raw = fopen(curr_data_raw_filename, 'r');

    % Loop through data chunks and apply common average referencing
    chunkInd = 1;
    data_eof = feof(curr_data_raw_filename);
    while ~data_eof

        % Load in current data chunk
        curr_data = fread(fid_raw, [384 chunkSize], '*int16');

        % Subtract median across channels
        curr_data_centered = curr_data - median(curr_data,2);

        % Loop through ADC channel indicies
        curr_data_car = zeros(size(curr_data),'int16');
        for curr_adc_idx = 1:size(ch_adc_idx,1)
            % Get same-index channels across ADCs
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
        end

        fwrite(fid_car, curr_data_car, 'int16');

        ap.print_progress_fraction(chunkInd,nChunksTotal);
        
        chunkInd = chunkInd+1;
        data_eof = feof(curr_data_raw_filename);
        
    end

    % Close raw file
    fclose(fid_raw);

end

% Close CAR file
fclose(fid_car);

disp('Finished CAR');





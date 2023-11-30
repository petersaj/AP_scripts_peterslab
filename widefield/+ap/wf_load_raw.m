function im = wf_load_raw(animal,rec_day,rec_time,frames)
% wf_load_raw(data_dir,frames)
%
% Load raw frames from binary (captured with plab.rig.widefield)

data_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
if ~exist(data_path,'dir')
    error('Data path not found: %s',data_path);
end

%% Get metadata

% Get filename
metadata_fn = fullfile(data_path,sprintf('widefield_%s_metadata.bin',rec_time));

% Load metadata
% Metadata format: 9 x n frames
% [image height, image width, frame number, timestamp (y,m,d,H,M,S)]
curr_metadata_fid = fopen(metadata_fn,'r');
curr_metadata = reshape(fread(curr_metadata_fid,'double'),9,[]);
fclose(curr_metadata_fid);

% Get image size
im_size = reshape(curr_metadata(1:2,1),1,[]);

%% Load images

% Open file for reading
data_fn = fullfile(data_path,sprintf('widefield_%s_data.bin',rec_time));
data_fid = fopen(data_fn,'r');

% Load all selected frames
im = zeros([im_size,length(frames)],'uint16');
for curr_frame_idx = 1:length(frames)
    curr_frame_location = prod(im_size)*(frames(curr_frame_idx)-1)*2; % uint16: *2 bytes
    fseek(data_fid,curr_frame_location,-1);
    im(:,:,curr_frame_idx) = reshape(fread(data_fid,prod(im_size),'uint16=>uint16'),im_size);
end








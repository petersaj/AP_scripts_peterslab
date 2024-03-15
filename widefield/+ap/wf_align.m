function [im_aligned,im_tform] = wf_align(im_unaligned,animal,day,align_type,master_align)
% im_aligned = align_widefield(im_unaligned,animal,day,align_type,master_align)
%
% Align widefield images across days and animals
%
% im_unaligned - unaligned images
% animal - animal corresponding to im_unaligned, e.g. 'XX001'
% day - 'yyyy-mm-dd' corresponding to im_unaligned (cell array if multiple)
% align_type - type of alignment:
% master_align - used with 'new_animal': master template for alignment
% > default (otherwise/empty) - apply saved day/animal transform to image(s)
% > 'day_only' - apply only day (only) transform to image(s)
% > 'animal_only' - apply only animal transform (e.g. already day-aligned)
% > 'new_days' - make new alignment within animal across days
% > 'new_animal' - make new alignment from animal to master
% > 'create_master' - create new master alignment reference (ONCE!!)
% > 'create_submaster' - create sub-master (align to master)
% master_align - if 'new_animal', use this instead of master_vfs
%
% USAGE:
% 1) Create master alignment reference ('create_master') from average aligned VFSs from multiple animals
% 2) Across-day align ('new_days') using vasulature from all days of one animal
% 3) Across-animal align ('new_animal) using the aligned average VFS from new animal
% 4) Apply saved alignments to any images (specifying animal and day source)


%% Initialize

% Set filename for saved alignments
alignment_path = fullfile(plab.locations.server_path,'Users','Andy_Peters','widefield_alignment');
alignment_path_animal = fullfile(alignment_path,'animal_alignment');
alignment_filename = fullfile(alignment_path_animal,sprintf('wf_alignment_%s.mat',animal));

if ~exist(alignment_path,'dir')
    mkdir(alignment_path)
end
if ~exist(alignment_path_animal,'dir')
    mkdir(alignment_path_animal)
end

% Set empty align_type if unassigned
if ~exist('align_type','var')
    align_type = '';
end

switch align_type

    case 'new_days'
        %% Create across-day alignments for animal
        % (only run once for each animal - aligns to iterative average)

        % If no unaligned images, find and load all average images
        if isempty(im_unaligned)
            recordings = plab.find_recordings(animal);
            wf_days_idx = find(cellfun(@(x) any(x),{recordings.widefield})');

            day = {recordings(wf_days_idx).day};
            im_unaligned = cell(size(wf_days_idx));
            for curr_day_idx = 1:length(wf_days_idx)
                curr_day = recordings(wf_days_idx(curr_day_idx)).day;
                curr_avg_im_fn = plab.locations.filename('server', ...
                    animal,curr_day,[],'widefield','meanImage_blue.npy');
                im_unaligned{curr_day_idx} = single(readNPY(curr_avg_im_fn));
            end
        end

        % Set output size as the largest image
        [im_y,im_x] = cellfun(@size,im_unaligned);

        im_y_max = max(im_y);
        im_x_max = max(im_x);
        ref_size = [im_y_max,im_x_max];

        im_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(im) padarray(im, ...
            [im_y_max - size(im,1),im_x_max - size(im,2)]./2,0,'both'), ...
            im_unaligned,'uni',false),1,1,[]));

        % Intensity-normalize the image used for alignment
        im_unaligned_pad_norm = imgaussfilt(im_unaligned_pad,2)./ ...
            imgaussfilt(im_unaligned_pad,10);

        %         (set transform optimizer)
        % (for OnePlusOneEvolutionary)
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-7;
        optimizer.InitialRadius = 1e-5;
        optimizer.Epsilon = 1.5e-5;

        %         % (for RegularStepGradientDescent)
        %         [optimizer, metric] = imregconfig('monomodal');
        %         optimizer.GradientMagnitudeTolerance = 1e-7;
        %         optimizer.MaximumIterations = 300;
        %         optimizer.MaximumStepLength = 1e-5;
        %         optimizer.MinimumStepLength = 1e-7;
        %         optimizer.RelaxationFactor = 0.6;

        % (first pass: rigid transform to day 1)
        % (increasing the PyramidLevels made the biggest improvement)
        disp('Rigid aligning images...')
        im_ref = im_unaligned{1};
        im_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
        tform_matrix = cell(size(im_unaligned));
        for curr_im = 1:length(im_unaligned)
            im_tform = imregtform(im_unaligned{curr_im}, ...
                im_ref,'rigid',optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(im_unaligned{curr_im}, ...
                im_tform,'Outputview',imref2d(ref_size));
            tform_matrix{curr_im} = im_tform.T;
            im_aligned(:,:,curr_im) = curr_im_reg;
        end

        % Check and manually align bad alignments by control points
        while true
            f = AP_imscroll(im_aligned);
            clim(prctile(im_aligned(:),[1,99]));
            axis image;
            user_input = input('Manual align? (#,s=save, q=quit): ','s');
            if ishandle(f)
                close(f)
            end

            % Convert user input to number if numeric
            if ~isnan(str2double(user_input))
                user_input = str2double(user_input);
            end

            if ischar(user_input) && ismember(user_input,{'s','q'})
                % User selected save or quit
                break

            elseif isnumeric(user_input) && ismember(user_input,1:size(im_aligned,3))
                % User selected a day to manually correct

                % Select control points for selected image
                [movingPoints,fixedPoints] = cpselect( ...
                    mat2gray(im_aligned(:,:,user_input),prctile(reshape(im_aligned(:,:,user_input),[],1),[0,98])), ...
                    mat2gray(im_ref,double(prctile(im_ref(:),[0,98]))),'Wait',true);
                cp_tform = fitgeotform2d(movingPoints,fixedPoints,'similarity');
                cp_rigid_tform = rigidtform2d(cp_tform.RotationAngle,cp_tform.Translation);

                % Update image with new transform
                im_aligned(:,:,user_input) = ...
                    imwarp(im_aligned(:,:,user_input), ...
                    cp_rigid_tform,'Outputview',imref2d(ref_size));

                % Update combined transform matrix for image
                tform_matrix{user_input} = tform_matrix{user_input}*cp_rigid_tform.T;

            elseif isnumeric(user_input) && ~ismember(user_input,1:size(im_aligned,3))
                % User selected an invalid day number
                fprintf('Invalid day number selected: %d',user_input);
                continue

            else
                % User selected an invalid character
                fprintf('Invalid input: %s',user_input);
                continue

            end
        end

        % If user selected save, save transform matrix into structure
        if strcmp(user_input,'s')
            wf_tform = struct;
            wf_tform.animal = animal;
            wf_tform.day = reshape(day,[],1);
            wf_tform.day_tform = tform_matrix;
            wf_tform.ref_size = ref_size;

            save(alignment_filename,'wf_tform');
            disp(['Saved day transforms for ' animal '.'])
        elseif strcmp(user_input,'q')
            disp('Alignment not saved');
        end

    case 'new_animal'
        %% Create across-animal alignment for animal to master VFS

        if ~exist('master_align') || isempty(master_align)
            % Align master VFS by default
            disp('Aligning animal VFS to master VFS...')
            % Load master VFS
            master_vfs_fn = fullfile(alignment_path,'master_vfs.mat');
            load(master_vfs_fn);
            master_align = master_vfs;
        else
            disp('Aligning animal image to master image')
        end

        % Align animal image to master image
        ref_size = size(master_align);

        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-3;

        im_tform = imregtform(im_unaligned,master_align,'similarity',optimizer,metric);
        im_aligned = imwarp(im_unaligned,im_tform,'Outputview',imref2d(ref_size));

        figure;
        tiledlayout(2,2,'TileSpacing','tight');
        nexttile;
        imshowpair(master_align,im_unaligned);
        title('Unaligned');
        nexttile;
        imshowpair(master_align,im_aligned);
        title('Aligned');
        nexttile;
        imagesc(master_align);
        axis image off
        clim([-1,1]);
        colormap(AP_colormap('BWR'));
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
        title('Master')
        nexttile;
        imagesc(im_aligned);
        axis image off
        clim([-1,1]);
        colormap(AP_colormap('BWR'));
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
        title('Aligned')

        % Save transform matrix into previously saved alignment file
        if ~exist(alignment_filename,'file')
            error('%s: no day alignments saved',animal);
        end

        load(alignment_filename);

        confirm_save = strcmp(input( ...
            ['Save animal alignment for ' animal '? (y/n): '], ...
            's'),'y');
        if confirm_save
            wf_tform.animal_tform = im_tform.T;
            wf_tform.master_ref_size = ref_size;
            save(alignment_filename,'wf_tform');
            disp(['Saved new animal alignment for ' animal '.'])
        else
            disp('Not overwriting.')
        end

    case 'create_master'
        %% Create master VFS
        % Whenever this is run, any dependent analysis must be re-run
%         error('Really create new master? Manually override');

        disp('Creating new master VFS...')

        % Check inputs
        if ~iscell(im_unaligned)
            error('Expected cell array of unaligned VFS''s');
        end

        % Check with user: are L and R hemi same or opposite chirality?
        LR_symmetry = questdlg('Are left/right hemisphere colors same or opposite?','Choose VFS type','symmetric','opposite','symmetric');
        switch LR_symmetry
            case 'symmetric'
                side_sign = 1;
            case 'opposite'
                side_sign = -1;
        end
        % Pad and concatenate unaligned VFS's
        [vfs_y,vfs_x] = cellfun(@size,im_unaligned);

        vfs_y_max = max(vfs_y);
        vfs_x_max = max(vfs_x);

        vfs_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(vfs) padarray(vfs, ...
            [vfs_y_max - size(vfs,1),vfs_x_max - size(vfs,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));

        % Set alignment parameters
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;

        % Iterate averaging the VFS and the aligning
        disp('Iterating aligning average VFS...')

        vfs_ref = nanmean(vfs_unaligned_pad,3);
        ref_size = size(vfs_ref);

        n_loops = 5;
        for curr_loop = 1:n_loops

            vfs_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
            for curr_animal = 1:length(im_unaligned)
                im_tform = imregtform(im_unaligned{curr_animal}, ...
                    vfs_ref,'affine',optimizer,metric,'PyramidLevels',5);
                curr_im_reg = imwarp(im_unaligned{curr_animal}, ...
                    im_tform,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = im_tform.T;
                vfs_aligned(:,:,curr_animal) = curr_im_reg;
            end

            vfs_ref = nanmean(vfs_aligned,3);
            AP_print_progress_fraction(curr_loop,n_loops);
        end

        % Symmetrize the average aligned VFS
        disp('Symmetrizing aligned average VFS...')

        n_loops = 15;
        ref_size = size(vfs_ref);

        ref_reg = nan(size(vfs_ref,1),size(vfs_ref,2),n_loops);
        for curr_loop = 1:n_loops

            ref_im_symm = (vfs_ref + side_sign*fliplr(vfs_ref))./2;

            im_tform = imregtform(vfs_ref,ref_im_symm,'rigid', ...
                optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(vfs_ref,im_tform,'Outputview',imref2d(ref_size));
            tform_matrix = im_tform.T;
            vfs_ref = curr_im_reg;

            ref_reg(:,:,curr_loop) = vfs_ref;
            AP_print_progress_fraction(curr_loop,n_loops);
        end

        % Set make symmetric-average map the master
        symm_vfs = ref_reg(:,:,end);
        symm_vfs_mirror_avg = (symm_vfs + side_sign*fliplr(symm_vfs))./2;

        master_vfs = symm_vfs_mirror_avg;

        figure;
        imagesc(master_vfs);
        axis image off
        colormap(ap.colormap('RWB'));
        clim([-1,1]);
        title('Master VFS');

        % Save the master
        confirm_save = input('Save new master VFS? (y/n): ','s');
        if strcmp(confirm_save,'y')
            [master_vfs_filename,master_vfs_path] = ...
                uiputfile([alignment_path filesep '*.mat'],'Save master VFS');
            master_vfs_fn = fullfile(master_vfs_path,master_vfs_filename);
            save(master_vfs_fn,'master_vfs');
            disp(['Saved new master VFS: ' master_vfs_fn]);
        else
            disp('Not saved.')
        end

    case 'create_submaster'
        %% Make submaster for new type, aligned to master VFS

        disp('Creating new submaster VFS...')

        % Check inputs
        if ~iscell(im_unaligned)
            error('Expected cell array of unaligned VFS''s');
        end

        % Check with user: are L and R hemi same or opposite chirality?
        LR_symmetry = questdlg('Are left/right hemisphere colors same or opposite?','Choose VFS type','symmetric','opposite','symmetric');
        switch LR_symmetry
            case 'symmetric'
                side_sign = 1;
            case 'opposite'
                side_sign = -1;
        end

        % Pad and concatenate unaligned VFS's
        [vfs_y,vfs_x] = cellfun(@size,im_unaligned);

        vfs_y_max = max(vfs_y);
        vfs_x_max = max(vfs_x);

        vfs_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(vfs) padarray(vfs, ...
            [vfs_y_max - size(vfs,1),vfs_x_max - size(vfs,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));

        % Set alignment parameters
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;

        % Iterate averaging the VFS and the aligning
        disp('Iterating aligning average VFS...')

        vfs_ref = nanmean(vfs_unaligned_pad,3);
        ref_size = size(vfs_ref);

        n_loops = 5;
        for curr_loop = 1:n_loops

            vfs_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
            for curr_animal = 1:length(im_unaligned)
                im_tform = imregtform(im_unaligned{curr_animal}, ...
                    vfs_ref,'affine',optimizer,metric,'PyramidLevels',5);
                curr_im_reg = imwarp(im_unaligned{curr_animal}, ...
                    im_tform,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = im_tform.T;
                vfs_aligned(:,:,curr_animal) = curr_im_reg;
            end

            vfs_ref = nanmean(vfs_aligned,3);
            AP_print_progress_fraction(curr_loop,n_loops);
        end

        % Symmetrize the average aligned VFS
        disp('Symmetrizing aligned average VFS...')

        n_loops = 15;
        ref_size = size(vfs_ref);

        ref_reg = nan(size(vfs_ref,1),size(vfs_ref,2),n_loops);
        for curr_loop = 1:n_loops

            ref_im_symm = (vfs_ref + side_sign*fliplr(vfs_ref))./2;

            im_tform = imregtform(vfs_ref,ref_im_symm,'rigid', ...
                optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(vfs_ref,im_tform,'Outputview',imref2d(ref_size));
            tform_matrix = im_tform.T;
            vfs_ref = curr_im_reg;

            ref_reg(:,:,curr_loop) = vfs_ref;
            AP_print_progress_fraction(curr_loop,n_loops);
        end

        % Set make symmetric-average map the master
        symm_vfs = ref_reg(:,:,end);
        symm_vfs_mirror_avg = (symm_vfs + side_sign*fliplr(symm_vfs))./2;

        % Load the master, align the submaster
        alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
        master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
        load(master_vfs_fn);
        % (if hemispheres opposite, flip left side)
        if strcmp(LR_symmetry,'opposite')
            master_vfs(:,1:round(size(master_vfs,2)/2)) = ...
                master_vfs(:,1:round(size(master_vfs,2)/2))*-1;
        end

        im_tform = imregtform(symm_vfs_mirror_avg, ...
            master_vfs,'affine',optimizer,metric,'PyramidLevels',5);
        submaster_vfs = imwarp(symm_vfs_mirror_avg, ...
            im_tform,'Outputview',imref2d(size(master_vfs)));

        % Plot submaster and alignment
        figure;
        subplot(1,3,1);
        imagesc(master_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
        title('Master');
        subplot(1,3,2);
        imagesc(submaster_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
        title('New aligned submaster');
        subplot(1,3,3);
        imshowpair(abs(master_vfs),abs(submaster_vfs));
        title('Master/submaster alignment');

        % Save the submaster
        confirm_save = input(['Save new submaster VFS? (y/n): '],'s');
        if strcmp(confirm_save,'y')
            [master_vfs_filename,master_vfs_path] = ...
                uiputfile([alignment_path filesep '*.mat'],'Save submaster VFS');
            master_vfs_fn = [master_vfs_path master_vfs_filename];
            save(master_vfs_fn,'submaster_vfs');
            disp(['Saved new submaster VFS: ' master_vfs_fn]);
        else
            disp('Not saved.')
        end


    otherwise
        %% Apply alignments to data

        % Load animal alignment
        if exist(alignment_filename,'file')
            load(alignment_filename)
        else
            error('No alignments found for %s', animal);
        end

        % Check for day alignment
        curr_day_idx = strcmp(day,wf_tform.day);
        if ~any(curr_day_idx) && ~strcmp(align_type,'animal_only')
            error('No %s alignment found for %s',animal,day);
        end

        switch align_type
            case 'day_only'
                % Align day only (ignore animal, even if present)
                curr_tform = wf_tform.day_tform{curr_day_idx};
                ref_size = wf_tform.ref_size;
            case 'animal_only'
                % Align animal only (used if already day-aligned)
                curr_tform = wf_tform.animal_tform;
                ref_size = wf_tform.master_ref_size;
            otherwise
                % Apply both day and animal alignments (if available)
                curr_day_tform = wf_tform.day_tform{curr_day_idx};
                if ~isempty(wf_tform.animal_tform)
                    curr_animal_tform = wf_tform.animal_tform;
                    curr_tform = curr_day_tform*curr_animal_tform;
                    ref_size = wf_tform.master_ref_size;
                else
                    warning([animal ' ' day ': No animal alignment']);
                    curr_tform = curr_day_tform;
                    ref_size = wf_tform.ref_size;
                end
        end

        % Transform image, cast to input data type
        tform = affinetform2d;
        tform.T = curr_tform;
        im_aligned = cast(imwarp(im_unaligned,tform,'Outputview',imref2d(ref_size)), ...
            class(im_unaligned));


end

% %% Align CCF to VFS (RUN ONCE)
% 
% alignment_path = fullfile(plab.locations.server_path,'Users','Andy_Peters','widefield_alignment');
% 
% % Load master VFS
% master_vfs_fn = fullfile(alignment_path,'master_vfs.mat');
% load(master_vfs_fn);
% master_align = master_vfs;
% 
% % Load CCF structure tree
% allen_atlas_path = fileparts(which('template_volume_10um.npy'));
% st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));
% 
% % Load dorsal CCF areas
% load(fullfile(alignment_path,'dorsal_cortex_borders.mat'));
% 
% % Color visual areas by VFS by Zhuang/Waters eLife 2017 Fig 3c
% dorsal_areas_unique = unique(dorsal_ccf_annotation(dorsal_ccf_annotation ~= 0));
% a_idx = find(cellfun(@(name) contains(name,'Anterior area'),st.safe_name(dorsal_areas_unique)));
% al_idx = find(cellfun(@(name) contains(name,'Anterolateral visual area'),st.safe_name(dorsal_areas_unique)));
% am_idx = find(cellfun(@(name) contains(name,'Anteromedial visual area'),st.safe_name(dorsal_areas_unique)));
% lm_idx = find(cellfun(@(name) contains(name,'Lateral visual area'),st.safe_name(dorsal_areas_unique)));
% v1_idx = find(cellfun(@(name) contains(name,'Primary visual area'),st.safe_name(dorsal_areas_unique)));
% p_idx = find(cellfun(@(name) contains(name,'Posterolateral visual area'),st.safe_name(dorsal_areas_unique)));
% pm_idx = find(cellfun(@(name) contains(name,'posteromedial visual area'),st.safe_name(dorsal_areas_unique)));
% li_idx = find(cellfun(@(name) contains(name,'Laterointermediate area'),st.safe_name(dorsal_areas_unique)));
% rl_idx = find(cellfun(@(name) contains(name,'Rostrolateral area'),st.safe_name(dorsal_areas_unique)));
% 
% ccf_vfs = zeros(size(dorsal_ccf_annotation));
% ccf_vfs(ismember(dorsal_ccf_annotation,dorsal_areas_unique( ...
%     [v1_idx,am_idx,al_idx,li_idx]))) = -1;
% ccf_vfs(ismember(dorsal_ccf_annotation,dorsal_areas_unique( ...
%     [a_idx,p_idx,pm_idx,rl_idx,lm_idx]))) = 1;
% 
% % Threshold VFS
% % (reference for alignment: Waters/Thompson, PLoS ONE 2019)
% vfs_cutoff = 0.05;
% master_vfs_thresh = zeros(size(master_vfs));
% master_vfs_thresh(master_vfs < -vfs_cutoff) = -1;
% master_vfs_thresh(master_vfs > vfs_cutoff) = 1;
% 
% % % (auto-align)
% % [optimizer, metric] = imregconfig('monomodal');
% % optimizer = registration.optimizer.OnePlusOneEvolutionary();
% % optimizer.MaximumIterations = 200;
% % optimizer.GrowthFactor = 1+1e-3;
% % optimizer.InitialRadius = 1e-6;
% % ccf_tform = imregtform(ccf_vfs,master_vfs_thresh,'affine',optimizer,metric,'PyramidLevels',5);
% % (control point align)
% [movingPoints,fixedPoints] = cpselect( ...
%     mat2gray(ccf_vfs),mat2gray(master_vfs_thresh),'Wait',true);
% ccf_tform = fitgeotform2d(movingPoints,fixedPoints,'affine');
% 
% dorsal_cortex_borders_aligned_long = cellfun(@(areas) cellfun(@(coords) ...
%     [fliplr(coords),ones(size(coords,1),1)]*ccf_tform.T,areas,'uni',false), ...
%     dorsal_cortex_borders,'uni',false);
% dorsal_cortex_borders_aligned = cellfun(@(areas) cellfun(@(coords) ...
%     coords,areas,'uni',false),dorsal_cortex_borders_aligned_long,'uni',false);
% 
% % Plot alignment
% figure;
% subplot(1,2,1);
% imagesc(ccf_vfs); hold on; axis image off
% cellfun(@(areas) cellfun(@(outline) ...
%     plot(outline(:,2),outline(:,1),'color','k'),areas,'uni',false), ...
%     dorsal_cortex_borders,'uni',false);
% colormap(AP_colormap('BWR'));
% 
% subplot(1,2,2);
% imagesc(master_vfs_thresh); hold on; axis image off
% cellfun(@(areas) cellfun(@(outline) ...
%     plot(outline(:,1),outline(:,2),'color','k'),areas,'uni',false), ...
%     dorsal_cortex_borders_aligned,'uni',false);
% colormap(AP_colormap('BWR'));
% 
% % Save master alignment
% confirm_save = input('Overwrite CCF to master VFS alignment? (y/n): ','s');
% if strcmp(confirm_save,'y')
%     save_ccf_fn = fullfile(alignment_path,'dorsal_cortex_borders_aligned');
%     save_ccf_tform_fn = fullfile(alignment_path,'ccf_tform.mat');
% 
%     save(save_ccf_fn,'dorsal_cortex_borders_aligned');
%     save(save_ccf_tform_fn,'ccf_tform');
%     disp('Saved CCF to master VFS alignment.')
% else
%     disp('Not saved.')
% end

















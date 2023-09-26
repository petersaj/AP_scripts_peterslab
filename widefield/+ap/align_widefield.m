function im_aligned = align_widefield(im_unaligned,animal,day,align_type,master_align)
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

% Set path with saved alignments
alignment_path = fullfile(plab.locations.server_path,'Users','Andy_Peters','widefield_alignment');
if ~exist(alignment_path,'dir')
    mkdir(alignment_path)
end

% Load transform structure
wf_tform_fn = fullfile(alignment_path,'wf_tform.mat');
% (if it doesn't exist yet, create it)
if ~exist(wf_tform_fn,'file')
    wf_tform = struct('animal',cell(0),'day',cell(0),'day_tform',cell(0), ...
        'ref_size',cell(0),'animal_tform',cell(0),'master_ref_size',cell(0));
    save(wf_tform_fn,'wf_tform');
else
    load(wf_tform_fn);
end

% Set empty align_type if unassigned
if ~exist('align_type','var')
    align_type = '';
end

%% Apply or create alignment

switch align_type
    
    case 'new_days'
        %% Align all days from one animal
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
                    animal,curr_day,[],'widefield','meanImage_violet.npy');
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
            tformEstimate_affine = imregtform(im_unaligned{curr_im}, ...
                im_ref,'rigid',optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(im_unaligned{curr_im}, ...
                tformEstimate_affine,'Outputview',imref2d(ref_size));
            tform_matrix{curr_im} = tformEstimate_affine.T;
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
                    mat2gray(im_aligned(:,:,user_input)), ...
                    mat2gray(im_ref),'Wait',true);
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
            curr_animal_idx = strcmp(animal,{wf_tform.animal});
            if isempty(curr_animal_idx) || ~any(curr_animal_idx)
                % (if not already extant, make new)
                curr_animal_idx = length(wf_tform) + 1;
                wf_tform(curr_animal_idx).animal = animal;
            end
            
            wf_tform(curr_animal_idx).day = reshape(day,[],1);
            wf_tform(curr_animal_idx).day_tform = tform_matrix;
            wf_tform(curr_animal_idx).ref_size = ref_size;
            
            save(wf_tform_fn,'wf_tform');
            disp(['Saved day transforms for ' animal '.'])
        elseif strcmp(user_input,'q')
            disp('Alignment not saved');
        end
        
        
    case 'create_master'
        %% Create master VFS
        % Whenever this is run, any dependent analysis must be re-run
        error('Really create new master? Manually override');
        
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
                tformEstimate_affine = imregtform(im_unaligned{curr_animal}, ...
                    vfs_ref,'affine',optimizer,metric,'PyramidLevels',5);
                curr_im_reg = imwarp(im_unaligned{curr_animal}, ...
                    tformEstimate_affine,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = tformEstimate_affine.T;
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
            
            tformEstimate_affine = imregtform(vfs_ref,ref_im_symm,'rigid', ...
                optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(vfs_ref,tformEstimate_affine,'Outputview',imref2d(ref_size));
            tform_matrix = tformEstimate_affine.T;
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
        colormap(brewermap([],'*RdBu'));
        caxis([-1,1]);
        title('Master VFS');
        
        % Save the master
        confirm_save = input(['Save new master VFS? (y/n): '],'s');
        if strcmp(confirm_save,'y')
            [master_vfs_filename,master_vfs_path] = ...
                uiputfile([alignment_path filesep '*.mat'],'Save master VFS');
            master_vfs_fn = [master_vfs_path master_vfs_filename];
            save(master_vfs_fn,'master_vfs');
            disp(['Saved new master VFS: ' master_vfs_fn]);
            confirm_align_ccf = strcmp(input(['Align CCF to new master VFS? (y/n): '],'s'),'y');
            if confirm_align_ccf
                AP_vfs_ccf_align;
            end
        else
            disp('Not saved.')
        end
        
    case 'create_submaster'
        %% Make submaster for new type, aligned to master
        
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
                tformEstimate_affine = imregtform(im_unaligned{curr_animal}, ...
                    vfs_ref,'affine',optimizer,metric,'PyramidLevels',5);
                curr_im_reg = imwarp(im_unaligned{curr_animal}, ...
                    tformEstimate_affine,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = tformEstimate_affine.T;
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
   
            tformEstimate_affine = imregtform(vfs_ref,ref_im_symm,'rigid', ...
                optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(vfs_ref,tformEstimate_affine,'Outputview',imref2d(ref_size));
            tform_matrix = tformEstimate_affine.T;
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
                
        tformEstimate_affine = imregtform(symm_vfs_mirror_avg, ...
            master_vfs,'affine',optimizer,metric,'PyramidLevels',5);
        submaster_vfs = imwarp(symm_vfs_mirror_avg, ...
            tformEstimate_affine,'Outputview',imref2d(size(master_vfs)));
        
        % Plot submaster and alignment
        figure; 
        subplot(1,3,1);
        imagesc(master_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title('Master');
        subplot(1,3,2);
        imagesc(submaster_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
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
        

    case 'new_animal'
        %% Align animal to master VFS
        
        if ~exist('master_align') || isempty(master_align)
            % Align master VFS by default          
            disp('Aligning animal VFS to master VFS...')
            % Load master VFS
            master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
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
        
        tformEstimate_affine = imregtform(im_unaligned,master_align,'similarity',optimizer,metric);
        im_aligned = imwarp(im_unaligned,tformEstimate_affine,'Outputview',imref2d(ref_size));
        tform_matrix = tformEstimate_affine.T;
        
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
        ap.draw_wf_ccf('ccf_aligned',[0.5,0.5,0.5]);
        title('Master')
        nexttile;
        imagesc(im_aligned);
        axis image off
        clim([-1,1]);
        colormap(AP_colormap('BWR'));
        ap.draw_wf_ccf('ccf_aligned',[0.5,0.5,0.5]);
        title('Aligned')
        
        % Save transform matrix into structure
        curr_animal_idx = strcmp(animal,{wf_tform.animal});
        if isempty(curr_animal_idx) || ~any(curr_animal_idx)
            % (if not already extant, make new)
            curr_animal_idx = length(wf_tform) + 1;
            confirm_save = true;
        else
            % (if extant, prompt to overwrite)
            confirm_save = strcmp(input(['Overwrite animal transform for ' animal '? (y/n): '],'s'),'y');
        end
        if confirm_save
            wf_tform(curr_animal_idx).animal = animal;
            wf_tform(curr_animal_idx).animal_tform = tform_matrix;
            wf_tform(curr_animal_idx).master_ref_size = ref_size;
            save(wf_tform_fn,'wf_tform');
            disp(['Saved new animal alignment for ' animal '.'])
        else
            disp('Not overwriting.')
        end
        

        
        
    otherwise
        %% Apply alignments to data
        
        % Find animal and day index within wf_tform structure
        curr_animal_idx = strcmp(animal,{wf_tform.animal});
        if ~any(curr_animal_idx)
           error(['No alignments found for ' animal]);
        end
        
        curr_day_idx = strcmp(day,wf_tform(curr_animal_idx).day);
        if ~any(curr_day_idx) && ~strcmp(align_type,'animal_only')
           error(['No ' animal ' alignment found for ' day]);
        end
             
        switch align_type
            case 'day_only'
                % Align day only (ignore animal, even if present)
                curr_tform = wf_tform(curr_animal_idx).day_tform{curr_day_idx};
                ref_size = wf_tform(curr_animal_idx).ref_size;
            case 'animal_only'
                % Align animal only (used if already day-aligned)
                curr_tform = wf_tform(curr_animal_idx).animal_tform;
                ref_size = wf_tform(curr_animal_idx).master_ref_size;
            otherwise
                % Apply both day and animal alignments (if available)
                curr_day_tform = wf_tform(curr_animal_idx).day_tform{curr_day_idx};
                if ~isempty(wf_tform(curr_animal_idx).animal_tform)
                    curr_animal_tform = wf_tform(curr_animal_idx).animal_tform;
                    curr_tform = curr_day_tform*curr_animal_tform;
                    ref_size = wf_tform(curr_animal_idx).master_ref_size;
                else
                    warning([animal ' ' day ': No animal alignment']);
                    curr_tform = curr_day_tform;
                    ref_size = wf_tform(curr_animal_idx).ref_size;
                end               
        end
            
        % Transform image, cast to input data type
        tform = affinetform2d;
        tform.T = curr_tform;
        im_aligned = cast(imwarp(im_unaligned,tform,'Outputview',imref2d(ref_size)), ...
            class(im_unaligned));
        

end












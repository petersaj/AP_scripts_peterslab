% Create master U's as a common basis set for combining widefield data

%% Get U's from all valid animals

% Find animals with across-animal alignments
wf_alignment_path = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment');
animal_alignment_path = fullfile(wf_alignment_path,'animal_alignment');


wf_alignment_dir = dir(fullfile(animal_alignment_path,'wf_alignment_*.mat'));
use_animal_idx = false(size(wf_alignment_dir));
for curr_animal = 1:length(wf_alignment_dir)
    curr_alignment = load(fullfile(animal_alignment_path, ...
        wf_alignment_dir(curr_animal).name));
    if isfield(curr_alignment.wf_tform,'animal_tform')
        use_animal_idx(curr_animal) = true;
    end
end

animals = extractBetween({wf_alignment_dir.name},'wf_alignment_','.mat');

% Get U's from first day of each animal
animal_U = cell(length(animals),1);
for curr_animal = 1:length(animals)

    preload_vars = who;

    animal = animals{curr_animal};

    recordings = plab.find_recordings(animal);
    use_recording = find(cellfun(@any,{recordings.widefield}),1);
    rec_day = recordings(use_recording).day;
    rec_time = recordings(use_recording).recording{1};

    load_parts = struct;
    load_parts.widefield = true;
    try
        ap.load_recording;
        if ~load_parts.widefield_align
            continue
        end
    catch me
        continue
    end

    animal_U{curr_animal} = wf_U;

    AP_print_progress_fraction(curr_animal,length(animals));
    clearvars('-except',preload_vars{:});
end



%% Create master U as SVD across aligned animal U's, save

[U_U,~,~] = svd(cell2mat(cellfun(@(x) reshape(x,[],size(x,3)), ...
    animal_U(~cellfun(@isempty,animal_U)),'uni',false)'),'econ');

n_components = 2000;

U_master = reshape(U_U(:,1:n_components),size(animal_U{1},1), ...
    size(animal_U{1},2),n_components);

save_fn = fullfile(wf_alignment_path,'U_master.mat');
save(save_fn,'U_master');







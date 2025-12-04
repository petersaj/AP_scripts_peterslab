%% Replace named node across ALL Bonsai workflows in entire repo
%
% !! WARNING !!
%
% This may break existing workflows
% Only use this for thoroughly tested changes, and notify everyone on use

% Find all bonsai workflows 
bonsai_dir = dir('C:\Github\PetersLab_rigging\bonsai_workflows\**\*.bonsai');
bonsai_fns = cellfun(@(folder,name) string(fullfile(folder,name)), ...
    {bonsai_workflows.folder},{bonsai_workflows.name});

% Load node to replace (previously saved into text file)
node_replace_fn = "C:\Users\petersa\Desktop\UpdateTrialNumber_new.txt";
node_replace = readlines(node_replace_fn);

% Loop through bonsai workflows, replace node if found

for curr_bonsai_fn = bonsai_fns

    % Read current bonsai file
    curr_bonsai_txt = readlines(curr_bonsai_fn);

    % Find the start of named node (which is created 1 line before the name)
    node_name = 'UpdateTrialNumber';
    node_start = find(contains(curr_bonsai_txt,node_name))-1;

    if isscalar(node_start)

        % Get expression open/closed and cumulative unclosed
        % (WorkflowOutputs open and don't close expressions? So ignore those)
        expression_opens = contains(curr_bonsai_txt,'<Expression') & ~contains(curr_bonsai_txt,'WorkflowOutput');
        expression_closes = contains(curr_bonsai_txt,'</Expression');
        expressions_unclosed = cumsum(expression_opens) - cumsum(expression_closes);

        % Find end of node by returning to same number of unclosed nodes
        node_end = find(expression_closes & ...
            expressions_unclosed == expressions_unclosed(node_start-1) & ...
            (1:length(curr_bonsai_txt))' > node_start,1);

        % Overwrite bonsai text with replaced node
        bonsai_txt_new = vertcat(curr_bonsai_txt(1:node_start-1),node_replace,curr_bonsai_txt(node_end+1:end));
        writelines(bonsai_txt_new,curr_bonsai_fn);

        fprintf('Replaced node: %s\n',curr_bonsai_fn);

    elseif length(node_start) > 1
        error('Multiple nodes found: %s\n',curr_bonsai_fn)

    elseif isempty(node_start)
        continue
    end

end
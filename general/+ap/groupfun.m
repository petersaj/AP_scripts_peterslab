function grouped_data = groupfun(group_function,data,varargin)
% grouped_data = groupfun(data,dim1 group, dim2 group...)
%
% Perform function on grouped data
% 
% data - ND matrix of data
% dimN group - one grouping variable  for each dimension (if empty, assumes
% no grouping)
%
%
% CURRENTLY: 
% - only does nanmean
% - does a for loop if one grouping dimension (quick), does fully indexed
% accumarray if multiple grouping dimensions (very slow), there's probably
% a smarter way to do this

% Check that there's a group variable for each dimension
groups = varargin;
if length(groups) ~= ndims(data)
    error('Number groups needs to match number of dimensions');
end

groups_dim = find(~cellfun(@isempty,groups));

% Get "full" groups (replace empty groups with full index)
empty_group_idx = cellfun(@isempty,groups);
groups_full = groups;
groups_full(empty_group_idx) = cellfun(@(x) (1:x)',num2cell(size(data,find(empty_group_idx))),'uni',false);

% Check that number of grouping variables matches dimension size
if ~isequal(cellfun(@length,groups_full),size(data))
    error('Group number doesn''t match dimension size')
end

% Get unique indicies for groups
% (sort from last column to first to keep column order after grouping)
[groups_unique,~,group_idx] = cellfun(@(x) unique(fliplr(x),'rows'),groups_full,'uni',false);

% Group data
if length(groups_dim) == 1
    % If only one grouping dimension:
    % concatenate along grouped dimension, loop through groups

    % Reshape data to ungrouped data x grouped data
    permute_order = [setdiff(1:ndims(data),groups_dim),groups_dim];
    [~,unpermute_order] = sort(permute_order);

    data = reshape(permute(data,permute_order),[],size(data,groups_dim));
    grouped_data = nan(size(data,1),max(group_idx{groups_dim}));

    for curr_grp = 1:max(group_idx{groups_dim})
        grouped_data(:,curr_grp) = nanmean(data(:,group_idx{groups_dim} == curr_grp),2);
    end

    % Reshape data to expected size and permutation
    grouped_data = permute(reshape(grouped_data,cellfun(@(x) size(x,1), ...
        groups_unique(permute_order))),unpermute_order);

else
    % If multiple grouping dimensions: 
    % create full grouping indicies and perform accumarray (slow)

    % Get full grid of group indicies
    group_grid = cell(size(groups_full));
    [group_grid{1:length(groups_full)}] = ndgrid(group_idx{:});

    % Perform accumarray with group indicies
    grouped_data = accumarray( ...
        cell2mat(cellfun(@(x) reshape(x,[],1),group_grid,'uni',false)), ...
        reshape(data,[],1), ...
        [], ...
        group_function, ...
        NaN);

end

% Reshape output to have each group on one dimension
grouped_data_size = cell2mat(cellfun(@(x) arrayfun(@(col) ...
    length(unique(x(:,col))),1:size(x,2)),groups_full,'uni',false));
grouped_data = reshape(grouped_data,grouped_data_size);















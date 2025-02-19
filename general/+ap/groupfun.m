function grouped_data = groupfun(group_function,data,varargin)
% grouped_data = groupfun(group_function, data, dim1 group, dim2 group...)
%
% Perform function on groups of data along any/multiple dimensions
% 
% INPUTS
% group_function: handle for function to perform on group
% data: ND matrix of data
% dimN group: grouping vector for each dimension, or empty if no group (if
% both data and group are vectors, can be any orientation). Groups that are
% NaN are ignored.
%
% OUTPUT
% grouped_data: ND data with function performed on group, with size reduced
% to number of groups in each grouped dimension.
%
% EXAMPLE: 
% Average values in groups along the 3rd and 4th dimensions
% rand_data = rand(3,4,5,6);
% group_3d = [1,1,1,2,2];
% group_4d = [1,1,2,2,3,3];
% grouped_data = ap.groupfun(@mean,rand_data,[],[],group_3d,group_4d)

% Check that there's a group variable for each dimension
groups = varargin;

if sum(size(data)~=1) ~= 1 && length(groups) ~= ndims(data)
    % If data has >1 non-singleton dimension, fill in empty groups for
    % number of dimensions (i.e. if no grouping specified for dimension,
    % don't group).
    groups{setdiff(1:ndims(data),find(~cellfun(@isempty,groups)))} = [];

elseif sum(size(data)~=1) == 1 && find(~cellfun(@isempty,groups)) == 1
    % If data has 1 non-singleton dimension and one group, put group in
    % non-singleton dimension slot
    data_orientation = find(size(data) ~= 1);
    groups{data_orientation} = groups{1};
    if length(groups) ~= 1
        groups(setdiff(1:length(groups),data_orientation)) = {[]};
    end

end

groups_dims = find(~cellfun(@isempty,groups));

% Check that number of grouping variables matches dimension size
if ~isequal(size(data,groups_dims), cellfun(@length, groups(groups_dims)))
    error('Group number doesn''t match dimension size')
end

% Permute data and reshaped to 2D: grouped x ungrouped
permute_order = [groups_dims,setdiff(1:ndims(data),groups_dims)];
[~,unpermute_order] = sort(permute_order);

data_reshape = reshape(permute(data,permute_order),prod(size(data,groups_dims)),[]);

% Get unique groups along combined group dimension
% (ignore entries with NaN group)
group_grid = cell(1,length(groups_dims));
[group_grid{:}] = ndgrid(groups{groups_dims});
group_flat = cell2mat(cellfun(@(x) reshape(x,[],1),group_grid,'uni',false));

% (flip groups so that the sorting hierarchy higher > lower)
group_flat_unique_idx = nan(size(group_flat,1),1);
group_nonan = ~any(isnan(group_flat),2);
[groups_flat_unique,~,group_flat_unique_idx(group_nonan)] = unique(fliplr(group_flat(group_nonan,:)),'rows');

% Loop through groups and do function
grouped_data = nan(max(group_flat_unique_idx),size(data_reshape,2),class(data));
for curr_grp = 1:max(group_flat_unique_idx)
    grouped_data(curr_grp,:) = group_function(data_reshape(group_flat_unique_idx == curr_grp,:),1);
end

% Reshape data to expected size and permutation
grouped_data_permuted_size = ...
    [cellfun(@(x) length(unique(x(~isnan(x)))),groups(groups_dims)), ...
    size(data,setdiff(1:ndims(data),groups_dims))];

grouped_data = permute(reshape(grouped_data,grouped_data_permuted_size), ...
    unpermute_order);











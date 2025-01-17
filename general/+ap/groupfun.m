function grouped_data = groupfun(group_function,data,varargin)
% grouped_data = groupfun(group_function, data, dim1 group, dim2 group...)
%
% Perform function on groups of data along any/multiple dimensions
% 
% INPUTS
% group_function: handle for function to perform on group
% data - ND matrix of data
% dimN group - grouping vector for each dimension, or empty if no group
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
if length(groups) ~= ndims(data)
    error('Number groups needs to match number of dimensions');
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
group_grid = cell(1,length(groups_dims));
[group_grid{:}] = ndgrid(groups{groups_dims});
group_flat = cell2mat(cellfun(@(x) reshape(x,[],1),group_grid,'uni',false));

% (flip groups so that the sorting hierarchy higher > lower)
[~,~,group_flat_unique] = unique(fliplr(group_flat),'rows');

% Loop through groups and do function
grouped_data = nan(max(group_flat_unique),size(data_reshape,2),class(data));
for curr_grp = 1:max(group_flat_unique)
    grouped_data(curr_grp,:) = group_function(data_reshape(group_flat_unique == curr_grp,:),1);
end

% Reshape data to expected size and permutation
grouped_data_permuted_size = ...
    [cellfun(@(x) length(unique(x)),groups(groups_dims)), ...
    size(data,setdiff(1:ndims(data),groups_dims))];

grouped_data = permute(reshape(grouped_data,grouped_data_permuted_size), ...
    unpermute_order);











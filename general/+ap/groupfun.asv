function [grouped_data,groups_output] = groupfun(group_function,data,varargin)
% grouped_data = groupfun(group_function, data, dim1 group, dim2 group...)
%
% Perform function on groups of data along any/multiple dimensions
% 
% INPUTS
% group_function: handle for function to perform on group
% data: ND matrix of data
% dimN group: grouping vector for each dimension, or empty if no group (if
% both data and group are vectors, can be any orientation). Groups that are
% NaN are ignored. Groups can have more than one label, in which case they
% are passed through 'unique' first. Multi-label groups must have labels in
% columns.
%
% OUTPUT
% grouped_data: ND data with function performed on group, with size reduced
% to number of groups in each grouped dimension.
% groups_output: groups corresponding to grouped data
%
% EXAMPLE: 
% Average values in groups along the 3rd and 4th dimensions
% rand_data = rand(3,4,5,6);
% group_3d = [1,1,1,2,2];
% group_4d = [1,1,2,2,3,3];
% grouped_data = ap.groupfun(@mean,rand_data,[],[],group_3d,group_4d)

% Check that there's a group variable for each dimension
groups = varargin;

if sum(size(data)>1) > 1 && length(groups) ~= ndims(data)
    % If data has >1 non-singleton dimension, fill in empty groups for
    % number of dimensions (i.e. if no grouping specified for dimension,
    % don't group).
    [groups{setdiff(1:ndims(data),find(~cellfun(@isempty,groups)))}] = deal([]);

elseif sum(size(data)~=1) == 1 && find(~cellfun(@isempty,groups)) == 1
    % If data has 1 non-singleton dimension and one group, put group in
    % non-singleton dimension slot
    data_orientation = find(size(data) ~= 1);
    groups{data_orientation} = groups{1};
    if length(groups) ~= 1
        groups(setdiff(1:length(groups),data_orientation)) = {[]};
    end

elseif isempty(data)
    % If no data entered, return empty outputs
    grouped_data = [];
    groups_output = [];
    return
end

groups_dims = find(~cellfun(@isempty,groups));

% Check that number of grouping variables matches dimension size
if ~isequal(size(data,groups_dims), cellfun(@length, groups(groups_dims)))
    error('Group number doesn''t match dimension size')
end

% If any groups are matricies, use as multi-label groups
% (labels must be in columns)
groups_sanitized = groups;
multilabel_group_idx = find(cellfun(@(x) sum(size(x)>1) > 1,groups));
multilabel_groups = cell(size(groups));
for curr_multilabel_group = multilabel_group_idx
    curr_group = groups{curr_multilabel_group};
    curr_group_unique = unique(curr_group(~any(isnan(curr_group),2),:),'rows','stable');
    [~,curr_group_renumbered] = ismember(curr_group,curr_group_unique,'rows');
    curr_group_renumbered(curr_group_renumbered == 0) = NaN;

    groups_sanitized{curr_multilabel_group} = curr_group_renumbered;
    multilabel_groups{curr_multilabel_group} = curr_group_unique;
end

% Permute data and reshaped to 2D: grouped x ungrouped
permute_order = [groups_dims,setdiff(1:ndims(data),groups_dims)];
[~,unpermute_order] = sort(permute_order);

data_reshape = reshape(permute(data,permute_order),prod(size(data,groups_dims)),[]);

% Get unique groups along combined group dimension
% (ignore entries with NaN group)
group_grid = cell(1,length(groups_dims));
[group_grid{:}] = ndgrid(groups_sanitized{groups_dims});
group_flat = cell2mat(cellfun(@(x) reshape(x,[],1),group_grid,'uni',false));

% (flip groups so that the sorting hierarchy higher > lower)
group_flat_unique_idx = nan(size(group_flat,1),1);
group_nonan = ~any(isnan(group_flat),2);
[groups_output_sanitized,~,group_flat_unique_idx(group_nonan)] = unique(fliplr(group_flat(group_nonan,:)),'rows');

% If there are multi-label groups: replace output groups with full group
if ~any(multilabel_group_idx)
    groups_output = groups_output_sanitized;
elseif any(multilabel_group_idx) && length(groups_dims) == 1
    groups_output = multilabel_groups{curr_multilabel_group}( ...
        groups_output_sanitized(:,groups_dims == curr_multilabel_group),:);
else
    error('Multiple groups with at least one multi-label group: output not organized yet');
    % Not sure whether to put full multi-label group in matrix or cell
end

% Loop through groups and do function
% (apply function to whole data if known dimension argument, otherwise
% apply function to data columns separately but slower)
grouped_data = nan(max(group_flat_unique_idx),size(data_reshape,2),class(data));
switch func2str(group_function)
    case {'mean','nanmean','median','nanmedian','sum','nansum'}
        % Functions with dimension as argument 2
        for curr_grp = 1:max(group_flat_unique_idx)
            grouped_data(curr_grp,:) = ...
                group_function(data_reshape(group_flat_unique_idx == curr_grp,:),1);
        end

    case {'std','nanstd','prctile','mad'}
        % Functions with dimension as argument 3
        for curr_grp = 1:max(group_flat_unique_idx)
            grouped_data(curr_grp,:) = ...
                group_function(data_reshape(group_flat_unique_idx == curr_grp,:),[],1);
        end

    otherwise
        % Function with no/unknown dimension argument: loop over columns
        for curr_grp = 1:max(group_flat_unique_idx)
            for curr_column = 1:size(data_reshape,2)
                grouped_data(curr_grp,curr_column) = ...
                    group_function(data_reshape(group_flat_unique_idx == curr_grp,curr_column));
            end
        end
end

% Reshape data to expected size and permutation
grouped_data_permuted_size = ...
    [cellfun(@(x) length(unique(x(~isnan(x)))),groups_sanitized(groups_dims)), ...
    size(data,setdiff(1:ndims(data),groups_dims))];

grouped_data = permute(reshape(grouped_data,grouped_data_permuted_size), ...
    unpermute_order);











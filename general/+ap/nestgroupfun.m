function [grouped_data,groups] = nestgroupfun(group_functions,data,group_labels,split_labels)
% [grouped_data,groups] = nestgroupfun(group_functions,data,group_labels,split_labels)
%
% Perform nested group functions
% NOTE: at the moment, this ONLY groups in the first dimension of 'data'
% 
% INPUTS 
% group_functions: cell array of figure handles for each level of grouping,
% N groups +1 (extra is for first-level grouping, which combines all data
% with the same labels across [group_labels,split_labels]).
% data: ND matrix of data to group
% group_labels: labels to group data in grouping order, size(data,1) x  N labels
% split labels: labels to retain data splits, size(data,1) x  N labels
%
% OUTPUTS
% grouped_data: data grouped hierarchically by 'group_labels'
% groups: group labels corresponding to grouped data, according to unique entires
% in 'split_labels'
%
% EXAMPLE Define groups by data halves, thirds, and odd/even. First,
% average all data with the same entries in all of these categories. Next,
% combine in order data by 1) averaging in the same half, then 2) median
% across odd or even. Retain the splitting of data by thirds.
% (define data and groups)
% x = (1:1000)';
% group_halves = max(1,ceil(linspace(0,2,length(x))'));
% group_thirds = max(1,ceil(linspace(0,3,length(x))'));
% group_odds = mod(x,2)+1;
% % (define grouping functions: number of grouping labels + 1)
% grouping_functions = {@mean,@mean,@median};
% % (define groups to combine, in nested order)
% group_labels = [group_halves,group_odds];
% % (define groups to retain split)
% split_labels = group_thirds;
% % (perform nested grouping)
% [grouped_data,groups] = ap.nestgroupfun(grouping_functions,x,group_labels,split_labels);


groups = [group_labels,split_labels];

groups_nest = cell(size(group_labels,2)+1,1);
groups_nest{1} = groups;

data_group_nest = cell(size(group_labels,2)+1,1);
data_group_nest{1} = data;

for curr_nest = 1:size(group_labels,2)+1

    curr_group = groups_nest{curr_nest};
    use_curr_group = ~any(isnan(curr_group),2);

    groups_nest_idx = nan(size(curr_group,1),1);
    [curr_group_nest,~,groups_nest_idx(use_curr_group)] = ...
        unique(curr_group(use_curr_group,:),'rows');

    if size(curr_group,1) > 1
        % If >1 data rows: group data
        data_group_nest{curr_nest+1} = ...
            ap.groupfun(group_functions{curr_nest},data_group_nest{curr_nest}, ...
            groups_nest_idx);
    elseif size(curr_group,1) == 1
        % If 1 data row: just carry forward row (can't group further)
        % (note: this is necessary because groupfun is flexible on
        % dimension arguments, so if a vector is input it assumes the
        % non-singleton dimension should be grouped)
        data_group_nest{curr_nest+1} =data_group_nest{curr_nest};
    end

    if curr_nest <= size(group_labels,2)
        % For all but last grouping: drop successive labels
        groups_nest{curr_nest + 1} = curr_group_nest(:,2:end);
    else
        % On last grouping: keep all group labels
        groups_nest{curr_nest + 1} = curr_group_nest;
    end

end

% Output final grouped data and groups
grouped_data = data_group_nest{end};
groups = groups_nest{end};











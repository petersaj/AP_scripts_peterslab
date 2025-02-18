function [grouped_data,groups] = nestgroupfun(group_functions,data,group_labels,split_labels)
%%% IN PROGRESS
% To finish in future: call groupfun in for loop across nested group labels

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

    data_group_nest{curr_nest+1} = ...
        ap.groupfun(group_functions{curr_nest},data_group_nest{curr_nest}, ...
        groups_nest_idx);

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











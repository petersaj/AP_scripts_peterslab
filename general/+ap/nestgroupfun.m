function grouped_data = groupfun(group_function,data,varargin)
% STARTED BUT DIDN'T FINISH: WASN'T ACTUALLY NEEDED FOR APPLICATION
%
% To finish in future: call groupfun in for loop across nested group labels


group_labels = [unit_animal];
split_labels = [unit_ld,unit_kidx];

groups = [group_labels,split_labels];
use_groups = ~any(isnan(groups),2);

groups_nest = cell(size(group_labels,2),1);
groups_nest{1} = groups;
for curr_nest = 1:size(group_labels,2)

    curr_group = groups_nest{curr_nest};
    use_curr_group = ~any(isnan(curr_group),2);

    groups_nest_idx = nan(size(curr_group,1),1);
    [curr_group_nest,~,groups_nest_idx(use_curr_group)] = ...
        unique(curr_group(use_curr_group,:),'rows');

    groups_nest{curr_nest + 1} = curr_group_nest(:,2:end);

end













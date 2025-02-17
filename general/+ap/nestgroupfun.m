function grouped_data = groupfun(group_function,data,varargin)
% grouped_data = groupfun(group_function, data, dim1 group, dim2 group...)
%
% Perform function on groups of data along any/multiple dimensions
% 
% INPUTS
% group_function: handle for function to perform on group
% data - ND matrix of data
% dimN group - grouping vector for each dimension, or empty if no group (if
% both data and group are vectors, can be any orientation)
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



groups = [unit_kidx,unit_ld,unit_animal];
use_groups = ~any(isnan(groups),2);

groups_nest_idx = nan(size(groups,1),1);
[groups_nest,~,groups_nest_idx(use_groups)] = unique(groups(use_groups,:),'rows');

unit_kidx_rec_grouped = ap.groupfun(@mean,unit_psth_cat{curr_stim},groups_nest_idx,[]);



% testing groupfun

x = rand(10,4,3);
grp1 = 1:10;
grp2 = [1,NaN,2,2];
grp3 = [1,2,2];

x2 = ap.groupfun(@mean,x,[],grp2,grp3);













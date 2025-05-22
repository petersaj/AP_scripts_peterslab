function binned_idx = quantile_bin(n_values,n_bins)
% Create index of N bins with equal element counts
% (somehow matlab doesn't have a function do this already?)
%
% e.g. to create 5 bins with equal elements across 13 total elements: 
% binned_idx = quantile_bin(13,5);
% binned_idx = [1,1,1,2,2,3,3,3,4,4,5,5,5] 
%
% n_values - number of values to create index for
% n_bins - number of bins

binned_idx = min(floor(linspace(1,n_bins+1,n_values))',n_bins);


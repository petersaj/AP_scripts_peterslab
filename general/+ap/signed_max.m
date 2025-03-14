function [max,idx] = signed_max(x,dim)
% [max,idx] = signed_max(x,dim)
%
% Signed maximum amplitude

% Get max of absolute value with linear and subscript index
[max_abs,idx_linear] = nanmax(abs(x),[],dim,'linear');
[~,idx] = nanmax(abs(x),[],dim);

% Apply sign of max amplitude
max = max_abs.*sign(x(idx_linear));





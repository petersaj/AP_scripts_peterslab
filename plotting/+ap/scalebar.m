function h = scalebar(varargin)
% ap.scalebar(x_scale,y_scale)
% ap.scalebar(axis,x_scale,y_scale)
% 
% Draws scalebars in the bottom left corner of a plot 

% First two numeric inputs are [x,y]
[x_scale,y_scale] = deal(varargin{cellfun(@isnumeric,varargin)});

% Set axis to draw on
if nargin == 3 && isgraphics(varargin{1})
    draw_ax = varargin{1};
else
    draw_ax = gca;
end

% Set properties
% (draw magenta for now so they're obvious when transferring)
scalebar_col = 'm';
scalebar_width = 2;

h = gobjects(2,1);
if ~isempty(x_scale)
    h(1) = line(draw_ax,min(xlim(draw_ax)) + [0,x_scale],repmat(min(ylim(draw_ax)),2,1),'color',scalebar_col,'linewidth',scalebar_width);
end
if ~isempty(y_scale)
    h(2) = line(draw_ax,repmat(min(xlim(draw_ax)),2,1),min(ylim(draw_ax)) + [0,y_scale],'color',scalebar_col,'linewidth',scalebar_width);
end


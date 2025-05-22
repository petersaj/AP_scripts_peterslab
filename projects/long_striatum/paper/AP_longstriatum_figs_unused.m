%% [Fig 1X] Cortex map depths pre/post learning

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

n_depths = 6;

ctx_map_depthmean = cellfun(@(x) ap.groupfun(@mean,x,[],[],discretize(1:size(x,3),n_depths)), ...
    all_ctx_maps_to_str.cortex_kernel_px,'uni',false);

ctx_map_depthmean_prelearn = nanmean(cat(4,ctx_map_depthmean{bhv.days_from_learning<0 & cellfun(@(x) size(x,3)==n_depths,ctx_map_depthmean)}),4);
ctx_map_depthmean_postlearn = nanmean(cat(4,ctx_map_depthmean{bhv.days_from_learning>=0 & cellfun(@(x) size(x,3)==n_depths,ctx_map_depthmean)}),4);
ap.imscroll(cat(4,ctx_map_depthmean_prelearn,ctx_map_depthmean_postlearn));
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));


figure;tiledlayout(n_depths,2,'TileSpacing','none');
colormap(ap.colormap('PWG',[],1.5));
for curr_depth = 1:n_depths
    nexttile;
    imagesc(ctx_map_depthmean_prelearn(:,:,curr_depth));
    axis image off;
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
    clim([-1,1]*0.05);

    nexttile;
    imagesc(ctx_map_depthmean_postlearn(:,:,curr_depth));
    axis image off;
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
    clim([-1,1]*0.05);
end
ap.prettyfig;
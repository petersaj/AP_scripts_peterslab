function [roi_trace,roi_mask] = wf_roi(U,V,guide_im,overlay,roi_mask)
% [roi_trace,roi_mask] = wf_roi(U,V,guide_im,overlay,roi_mask)
%
% Draw roi on average SVD'd image, return reconstructed trace in time
%
% U - SVD U (Y x X x components)
% V - SVD V (components x T x ...[conditions])
% guide_im - can be an image or 'master' if master-aligned CCF/retinotopy
% overlay - optional, if drawing ROI mask, overlays this image
% roi_mask - optional, used if provided and prompted to draw otherwise (can
% also be Y x X x N rois)

% Choose ROI
if ~exist('roi_mask','var') || isempty(roi_mask)
    h = figure;
    
    ax_guide_im = axes;
    if isnumeric(guide_im)
        imagesc(ax_guide_im,guide_im);
        axis image off;
        colormap(ax_guide_im, gray);
        caxis(ax_guide_im, [prctile(guide_im(:),0.5) prctile(guide_im(:),99.5)]);
    elseif ischar(guide_im) && strcmp(guide_im,'master')        
        imagesc(ax_guide_im,ones(size(U(:,:,1))));
        colormap(gray);
        axis image off;
        caxis([0,1])
        set(gca,'YDir','Reverse');        
        AP_reference_outline('retinotopy','r');
        AP_reference_outline('ccf_aligned','k');        
    end
    if exist('overlay','var') && ~isempty(overlay) && all(size(overlay) == size(guide_im))
        ax_overlay = axes;
        overlay_im = imagesc(ax_overlay,overlay);
        axis image off;
        colormap(ax_overlay, hsv);
        overlay_thresh = double(prctile(abs(overlay(:)),95));
        caxis(ax_overlay, [-overlay_thresh,overlay_thresh]);
        set(overlay_im,'AlphaData',0.3*mat2gray(abs(overlay),[0,overlay_thresh]));
    end
    
    title('Draw ROI');
    roi_mask = roipoly;
    close(h);
end


% Get fluorescence across session in ROIs
if exist('V','var') && ~isempty(V)
    V_size = size(V);
    
    % Weight and sum pixels in U
    U_roi = transpose(reshape(U,[],size(U,3))'*reshape(roi_mask,[],size(roi_mask,3)));            

    % Reconstruct ROI trace
    roi_trace = reshape(U_roi*reshape(V,size(V,1),[]),[size(U_roi,1),V_size(2:end)]);
    
    % If mask was binary, divide by n pixels to make trace ROI average
    % (if mask was weighted, then shouldn't divide)
    if all(ismember(roi_mask(:),[0,1]))    
        roi_trace = roi_trace./(sum(reshape(roi_mask,[],size(roi_mask,3),1))');
    end
    
else
    roi_trace = [];
end




function [x, y, raw] = digitize2D(imgPath, xlim_data, t_guess, N)
% DIGITIZE_AIRFOIL_XY
% Auto-digitize a red airfoil outline and return x,y (panel-ready boundary).
% 
% Inputs
%   imgPath   : image path (e.g., 'NACA0018.png')
%   xlim_data : [xmin xmax], usually [0 1]
%   t_guess   : thickness ratio (e.g., 0.18 for NACA 0018)
%   N         : (optional) number of boundary panels (default 1000)
%
% Outputs
%   x, y : ordered boundary (row vectors), TE(lower)->LE->TE(upper), length N+1
%   raw  : struct with raw scattered points (x_over_c, y_over_c)
%
% Notes
%   - Finds largest red region, traces its boundary (no thinning),
%     maps LE->x=0 and TE->x=1 from the curve itself,
%     recenters camber ~ 0 and scales thickness ~ t_guess,
%     then resamples onto cosine-spaced x for panel accuracy.

    if nargin < 4 || isempty(N), N = 1000; end

    % --- load image ---
    I0 = imread(imgPath);
    I  = im2double(I0);

    % --- segment RED curve (robust thresholds; tweak if needed) ---
    HSV = rgb2hsv(I);
    H = HSV(:,:,1); S = HSV(:,:,2); V = HSV(:,:,3);
    mask_red = ((H < 0.07) | (H > 0.93)) & (S > 0.25) & (V > 0.20);
    mask_red = bwareaopen(mask_red, 30);
    mask_red = imclose(mask_red, strel('disk',3));
    mask_red = imfill(mask_red, 'holes');

    % keep largest red blob
    CC = bwconncomp(mask_red);
    assert(CC.NumObjects>0,'No red region found – adjust thresholds or crop.');
    [~,bigIdx] = max(cellfun(@numel, CC.PixelIdxList));
    mask_big = false(size(mask_red));
    mask_big(CC.PixelIdxList{bigIdx}) = true;

    % --- trace boundary pixels (no thinning → preserves LE) ---
    B = bwboundaries(mask_big,'noholes');
    b = B{1};                           % [row=y, col=x]
    yPix = b(:,1);
    xPix = b(:,2);

    % --- map X using the red curve: LE->0, TE->1 ---
    LE_pix = min(xPix);                  % leftmost pixel
    TE_pix = max(xPix);                  % rightmost pixel
    assert(TE_pix>LE_pix,'Failed to infer LE/TE from boundary.');

    sx = (xlim_data(2) - xlim_data(1)) / (TE_pix - LE_pix);
    bx = xlim_data(1) - sx*LE_pix;

    % provisional x/c for binning
    x_over_c_tmp = sx*xPix + bx;

    % --- build upper/lower pixel envelopes vs x to center & scale thickness ---
    xb_fine = linspace(0,1,1001).';
    yu_pix = nan(size(xb_fine));  yl_pix = nan(size(xb_fine));
    dx = 0.001;                               % half-width in x/c for binning
    for k = 1:numel(xb_fine)
        in = abs(x_over_c_tmp - xb_fine(k)) < dx;
        if any(in)
            % image rows increase downward → upper has smaller y
            yu_pix(k) = min(yPix(in));
            yl_pix(k) = max(yPix(in));
        end
    end
    yu_pix = fillmissing(yu_pix,'linear');
    yl_pix = fillmissing(yl_pix,'linear');

    % center camber ≈ 0 and scale so thickness ≈ t_guess
    yc_pix = 0.5*(yu_pix + yl_pix);
    yu0 = yu_pix - yc_pix;
    yl0 = yl_pix - yc_pix;
    t_pix = mean(yl0 - yu0, 'omitnan');       % pixel thickness
    assert(t_pix>0,'Thickness detection failed.');
    sy = t_guess / t_pix;                     % data units per pixel

    % map boundary pixels to (x/c, y/c) — subtract local camber then scale
    yc_at = interp1(xb_fine, yc_pix, min(max(x_over_c_tmp,0),1), 'linear','extrap');
    x_raw = sx*xPix + bx;
    y_raw = sy * (-(yPix - yc_at));           % minus: pixel rows increase downward

    % keep in [0,1]
    kkeep = (x_raw>=xlim_data(1)) & (x_raw<=xlim_data(2));
    x_raw = x_raw(kkeep);
    y_raw = y_raw(kkeep);

    % return raw scattered too (optional)
    raw.x_over_c = x_raw(:);
    raw.y_over_c = y_raw(:);

    % --- make smooth upper/lower on a clean x-grid, recenter camber again ---
    xb = linspace(0,1,1001).';
    yu = nan(size(xb)); yl = nan(size(xb));
    for k = 1:numel(xb)
        in = abs(x_raw - xb(k)) < dx;
        if any(in)
            yu(k) = max(y_raw(in));          % upper y/c
            yl(k) = min(y_raw(in));          % lower y/c
        end
    end
    yu = fillmissing(yu,'linear'); 
    yl = fillmissing(yl,'linear');
    yc = 0.5*(yu+yl); yu = yu - yc; yl = yl - yc;

    % --- resample to cosine-spaced x and build panel-ready boundary ---
    beta = linspace(0,pi,N+1);
    x_cos = 0.5*(1 - cos(beta));             % 0..1
    yU = interp1(xb, yu, x_cos, 'pchip','extrap');
    yL = interp1(xb, yl, x_cos, 'pchip','extrap');

    % TE(lower) -> LE -> TE(upper)
    x_lower = fliplr(x_cos);  y_lower = fliplr(yL);
    x_upper = x_cos(2:end);   y_upper = yU(2:end);

    x = [x_lower, x_upper];                   % row vectors, length N+1
    y = [y_lower, y_upper];
end

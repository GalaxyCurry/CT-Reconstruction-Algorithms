function img=cropCBCT(img,geo,enableConeCrop)
% cropCBCT Crops CBCT reconstruction to the effective field of view
%
% img = cropCBCT(img, geo)
%     Applies only cylindrical cropping (default).
%
% img = cropCBCT(img, geo, true)
%     Applies cylindrical cropping + conical trimming on top and bottom.
%
% This function sets voxels outside the effective FOV to zero.
% The FOV is defined by the system geometry (DSO, DSD, detector size).

% Default: disable conical crop (only cylinder)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

if nargin < 3
    enableConeCrop = false;
end

% Tangent is equal, cropRadious:
cropRadiusPhysical=(geo.sDetector(1,1)/2*geo.DSO(1))/geo.DSD(1);
% maximum distance from O
maxD=min(geo.nVoxel(1,1))/2;

% Crop radious will be theminimum of them
cropRadius=min([cropRadiusPhysical/geo.dVoxel(1,1) maxD]);

[x,y]=meshgrid(1:size(img,1),1:size(img,2));
inM=(x-size(img,1)/2).^2+(y-size(img,2)/2).^2<cropRadius^2;
%crop
img=bsxfun(@times,img,inM);


%% -------------------------------------------------------------------------
% Optional: Apply conical trimming on top and bottom slices
%% -------------------------------------------------------------------------
if enableConeCrop
    % Conical crop parameters
    cropH=(geo.sDetector(2,1)/2*geo.DSO(1))/geo.DSD(1);
    maxZ=geo.sVoxel(3,1)/2;
    cropH=maxZ-cropH;
    cropH2=cropH*(geo.DSO(1)-geo.sVoxel(1,1)/2)/geo.DSO(1);
    
    % Z physical positions (centered at 0)
    zPositions = ((1:size(img,3)) - (size(img,3)+1)/2) * geo.dVoxel(3,1); 
    
    % Precompute squared distance from center (XY grid)
    [xx, yy] = meshgrid((1:size(img,1)) - (size(img,1)+1)/2, ...
                        (1:size(img,2)) - (size(img,2)+1)/2);
    distSq = xx.^2 + yy.^2;
    
    % Apply conical mask slice by slice
    for zIdx = 1:nZ
        zAbs = abs(zPositions(zIdx));

        if zAbs <= cropH2
            % Central region: full radius, no change
            continue;
        end

        % Reduce radius linearly toward top/bottom
        excessRatio = (zAbs - cropH2) / (maxZ - cropH2);
        currentR = max(0, cropRadius * (1 - excessRatio));

        % Apply mask to current slice
        sliceMask = distSq < currentR^2;
        img(:,:,zIdx) = img(:,:,zIdx) .* single(sliceMask);
    end
end



end
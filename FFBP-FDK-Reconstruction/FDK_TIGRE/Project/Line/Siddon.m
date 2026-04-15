function proj = Siddon(img, geo, angles)
% SIDDON
% Implementation of Siddon's ray-tracing algorithm for CBCT forward projection.
% 非等步长
% Inputs:
%   img     - 3D voxel volume [X, Y, Z]
%   geo     - TIGRE geometry structure
%   angles  - 3xN array of projection angles [alpha; theta; psi]
%
% Output:
%   proj    - Sinogram/projections [detectorV, detectorU, nAngles]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    nU      = geo.nDetector(1);
    nV      = geo.nDetector(2);
    nAng    = size(angles, 2);
    proj    = zeros(nV, nU, nAng, 'single');

    for iAng = 1:nAng
       
        fprintf('Processing projection %4d/%4d\n', iAng, nAng)

        % Compute source and detector geometry in voxel coordinates
        [source, detOrigin, du, dv] = computeGeometry(geo, iAng, angles);
        
        % Loop over detector pixels
        for v = 1:nV
            for u = 1:nU
                % Compute 3D world coordinate of current detector pixel
                x = detOrigin(1) + (u-1)*du(1) + (v-1)*dv(1);
                y = detOrigin(2) + (u-1)*du(2) + (v-1)*dv(2);
                z = detOrigin(3) + (u-1)*du(3) + (v-1)*dv(3);
                
                % Perform Siddon ray tracing
                proj(v, u, iAng) = siddon_core(img, geo, source, [x; y; z]);
            end
        end
    end
end

%% ========================================================================
function val = siddon_core(img, geo, src, pix)
% Core Siddon ray-voxel intersection algorithm
% Traces a ray from source to detector pixel and accumulates line integrals
%
% src   - Source position in voxel grid coordinates
% pix   - Detector pixel position in voxel grid coordinates

    nx = geo.nVoxel(1,1);
    ny = geo.nVoxel(2,1);
    nz = geo.nVoxel(3,1);

    % Ray direction vector
    rx = pix(1) - src(1);
    ry = pix(2) - src(2);
    rz = pix(3) - src(3);

    % Avoid division by zero
    eps = 1e-3;
    if abs(rx) < eps, rx = 0; end
    if abs(ry) < eps, ry = 0; end
    if abs(rz) < eps, rz = 0; end

    % Compute parametric entry/exit points with volume boundaries
    tx1 = (-src(1)) / rx;        tx2 = (nx - src(1)) / rx;
    ty1 = (-src(2)) / ry;        ty2 = (ny - src(2)) / ry;
    tz1 = (-src(3)) / rz;        tz2 = (nz - src(3)) / rz;

    t_enter = max( [min(tx1,tx2), min(ty1,ty2), min(tz1,tz2)] );
    t_exit  = min( [max(tx1,tx2), max(ty1,ty2), max(tz1,tz2)] );

    % No intersection with volume
    if t_enter >= t_exit
        val = 0;
        return;
    end

    % Starting voxel indices
    x0 = src(1) + t_enter * rx;
    y0 = src(2) + t_enter * ry;
    z0 = src(3) + t_enter * rz;

    i = floor(x0);
    j = floor(y0);
    k = floor(z0);

    % Step direction along each axis
    stepX = sign(rx);
    stepY = sign(ry);
    stepZ = sign(rz);

    % Parameters to next voxel boundary
    tX = nextT(src(1), rx, i, stepX);
    tY = nextT(src(2), ry, j, stepY);
    tZ = nextT(src(3), rz, k, stepZ);

    t = t_enter;
    total = 0.0;

    % Traverse voxels along the ray
    while t < t_exit
        t_next = min([tX, tY, tZ]);
        dt = t_next - t;
        
        % Accumulate contribution if inside valid volume
        if i >= 1 && i <= nx && j >= 1 && j <= ny && k >= 1 && k <= nz
            total = total + dt * img(i, j, k);
        end
        
        t = t_next;

        % Step to next voxel
        if tX == t_next
            i = i + stepX;
            tX = nextT(src(1), rx, i, stepX);
        elseif tY == t_next
            j = j + stepY;
            tY = nextT(src(2), ry, j, stepY);
        else
            k = k + stepZ;
            tZ = nextT(src(3), rz, k, stepZ);
        end
    end

    val = total;
end

%% ========================================================================
function t = nextT(s, r, idx, step)
% Compute parametric value t at which the ray crosses the next voxel boundary
% s     - Source coordinate along axis
% r     - Ray direction component
% idx   - Current voxel index
% step  - Direction of traversal (+1 or -1)

    if r == 0
        t = inf;
    else
        t = (idx + (step > 0) - s) / r;
    end
end

%% ========================================================================
function [src, det0, du, dv] = computeGeometry(geo, i, angles)
% Compute source and detector positions in normalized voxel coordinates
% Applies rotations, offsets, COR, scaling, and detector geometry

    DSO     = geo.DSO(1,i);
    DSD     = geo.DSD(1,i);
    dU      = geo.dDetector(1);
    dV      = geo.dDetector(2);
    nU      = geo.nDetector(1);
    nV      = geo.nDetector(2);

    offOrigX = geo.offOrigin(1,i);
    offOrigY = geo.offOrigin(2,i);
    offOrigZ = geo.offOrigin(3,i);
    offDetU  = geo.offDetector(1,i);
    offDetV  = geo.offDetector(2,i);
    cor      = geo.COR(1,i);

    dRoll    = geo.rotDetector(1,i);
    dPitch   = geo.rotDetector(2,i);
    dYaw     = geo.rotDetector(3,i);

    % Initial source position
    src = [DSO; 0; 0];

    % Initial detector corner and pixel steps
    det0 = [ -(DSD - DSO);
             dU*(-nU/2 + 0.5);
             dV*( nV/2 - 0.5) ];
    du = [0; dU; 0];
    dv = [0; 0; -dV];

    % Apply detector roll-pitch-yaw rotation
    det0 = rotRPY(det0, dRoll, dPitch, dYaw);
    du   = rotRPY(du,   dRoll, dPitch, dYaw);
    dv   = rotRPY(dv,   dRoll, dPitch, dYaw);

    % Apply detector offset
    det0(2) = det0(2) + offDetU;
    det0(3) = det0(3) + offDetV;

    % Apply global ZYZ Euler rotation
    det0 = rotZYZ(det0, angles(1,i), angles(2,i), angles(3,i));
    du   = rotZYZ(du,   angles(1,i), angles(2,i), angles(3,i));
    dv   = rotZYZ(dv,   angles(1,i), angles(2,i), angles(3,i));
    src  = rotZYZ(src,  angles(1,i), angles(2,i), angles(3,i));

    % Apply object (isocenter) offset
    det0(1) = det0(1) - offOrigX;
    det0(2) = det0(2) - offOrigY;
    det0(3) = det0(3) - offOrigZ;
    src(1)  = src(1)  - offOrigX;
    src(2)  = src(2)  - offOrigY;
    src(3)  = src(3)  - offOrigZ;

    % Translate to align volume origin to corner
    det0 = det0 + [geo.sVoxel(1,1)/2 - geo.dVoxel(1,1)/2;
                   geo.sVoxel(2,1)/2 - geo.dVoxel(2,1)/2;
                   geo.sVoxel(3,1)/2 - geo.dVoxel(3,1)/2];
    src = src + [geo.sVoxel(1,1)/2 - geo.dVoxel(1,1)/2;
                 geo.sVoxel(2,1)/2 - geo.dVoxel(2,1)/2;
                 geo.sVoxel(3,1)/2 - geo.dVoxel(3,1)/2];

    % Normalize coordinates to voxel units
    det0(1) = det0(1) / geo.dVoxel(1,1);
    det0(2) = det0(2) / geo.dVoxel(2,1);
    det0(3) = det0(3) / geo.dVoxel(3,1);

    du(1)   = du(1) / geo.dVoxel(1,1);
    du(2)   = du(2) / geo.dVoxel(2,1);
    du(3)   = du(3) / geo.dVoxel(3,1);

    dv(1)   = dv(1) / geo.dVoxel(1,1);
    dv(2)   = dv(2) / geo.dVoxel(2,1);
    dv(3)   = dv(3) / geo.dVoxel(3,1);

    src(1)  = src(1)  / geo.dVoxel(1,1);
    src(2)  = src(2)  / geo.dVoxel(2,1);
    src(3)  = src(3)  / geo.dVoxel(3,1);

    % Apply Center of Rotation (COR) correction
    cor_x = -cor * sin(angles(1,i)) / geo.dVoxel(1,1);
    cor_y =  cor * cos(angles(1,i)) / geo.dVoxel(2,1);
    det0(1) = det0(1) + cor_x;
    det0(2) = det0(2) + cor_y;
    src(1)  = src(1)  + cor_x;
    src(2)  = src(2)  + cor_y;
end

%% ========================================================================
function p = rotRPY(p, r, pch, y)
% Roll-Pitch-Yaw rotation (intrinsic Z-Y-X)
% Used for detector orientation

    R = [
        cos(r)*cos(pch),  cos(r)*sin(pch)*sin(y)-sin(r)*cos(y),  cos(r)*sin(pch)*cos(y)+sin(r)*sin(y);
        sin(r)*cos(pch),  sin(r)*sin(pch)*sin(y)+cos(r)*cos(y),  sin(r)*sin(pch)*cos(y)-cos(r)*sin(y);
        -sin(pch),        cos(pch)*sin(y),                       cos(pch)*cos(y)
    ];
    p = R * p;
end

%% ========================================================================
function p = rotZYZ(p, a, t, ps)
% ZYZ Euler angle rotation for gantry rotation

    R = [
        cos(a)*cos(t)*cos(ps)-sin(a)*sin(ps),  -cos(a)*cos(t)*sin(ps)-sin(a)*cos(ps),  cos(a)*sin(t);
        sin(a)*cos(t)*cos(ps)+cos(a)*sin(ps),  -sin(a)*cos(t)*sin(ps)+cos(a)*cos(ps),  sin(a)*sin(t);
        -sin(t)*cos(ps),                        sin(t)*sin(ps),                         cos(t)
    ];
    p = R * p;
end
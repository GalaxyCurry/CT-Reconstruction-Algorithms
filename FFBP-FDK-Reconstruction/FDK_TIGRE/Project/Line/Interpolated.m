function proj = Interpolated(img, geo, angles)
% TIGRE-style forward projection using 3D trilinear interpolation
% belongs to Joseph 投影 / 等步长采样 一类
%
% Inputs:
%   img     - 3D voxel array [nx, ny, nz]
%   geo     - CT geometry structure (TIGRE-compatible)
%   angles  - 3 x nAngles matrix [alpha, theta, psi]'
%
% Output:
%   proj    - Projection images [nDetecV, nDetecU, nAngles]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    
    nU   = geo.nDetector(1);
    nV   = geo.nDetector(2);
    nAng = size(angles, 2);
    proj = zeros(nV, nU, nAng, 'single');


    for iAng = 1:nAng

        fprintf('Processing angle %4d/%4d\n', iAng, nAng);
        
        
        % Compute geometry parameters: source, detector origin, pixel steps
        [source, detOrigin, du, dv] = computeDeltas(geo, iAng, angles);
        
        % Loop over detector pixels (u, v)
        for v = 1:nV  
            for u = 1:nU
                
                % 3D coordinate of current detector pixel
                x = detOrigin(1) + (u-1) * du(1) + (v-1) * dv(1);
                y = detOrigin(2) + (u-1) * du(2) + (v-1) * dv(2);
                z = detOrigin(3) + (u-1) * du(3) + (v-1) * dv(3);
                
                % Ray tracing along source -> detector pixel
                proj(v, u, iAng) = ray_trace_core(img, geo, source, [x; y; z], iAng);
            end
        end
    end
end

%% ============================================================================
function val = ray_trace_core(img, geo, source, P, iAng)
% Perform ray tracing with uniform step sampling + trilinear interpolation
% Integrate along the ray using fixed step size (geo.accuracy)
% Sum interpolated voxel values and multiply by physical step length

    % Ray direction vector
    dx = P(1) - source(1);
    dy = P(2) - source(2);
    dz = P(3) - source(3);

    ray_len = sqrt(dx^2 + dy^2 + dz^2);
    if ray_len < 1e-6
        val = 0;
        return;
    end

    % Total number of steps along the ray
    nSteps = ceil(ray_len / geo.accuracy);
    
    % Unit step in normalized voxel coordinates
    stepX = dx / nSteps;
    stepY = dy / nSteps;
    stepZ = dz / nSteps;

    % Skip safe distance before entering the volume 
    cropdist = floor(maxdistanceCuboid(geo, iAng) / geo.accuracy);
    sum_val  = 0.0;

    % Integrate along the ray
    for i = cropdist : nSteps
        x = source(1) + stepX * i;
        y = source(2) + stepY * i;
        z = source(3) + stepZ * i;

        % +0.5 offset matches CUDA texture sampling coordinate convention
        sum_val = sum_val + trilinear_interp(img, x+0.5, y+0.5, z+0.5);
    end

    % Physical length of one step in world coordinates
    dl = sqrt( (stepX * geo.dVoxel(1,1))^2 + ...
               (stepY * geo.dVoxel(2,1))^2 + ...
               (stepZ * geo.dVoxel(3,1))^2 );

    val = sum_val * dl;
end

%% ============================================================================
function val = trilinear_interp(img, x, y, z)

    [nx, ny, nz] = size(img);

    % Out-of-bounds → return zero (border mode in texture)
    if x < 1 || x > nx || y < 1 || y > ny || z < 1 || z > nz
        val = 0;
        return;
    end

    % Integer indices and fractional weights
    ix = floor(x); fx = x - ix;
    iy = floor(y); fy = y - iy;
    iz = floor(z); fz = z - iz;

    % Clamp to valid range
    ix = max(1, min(ix, nx-1));
    iy = max(1, min(iy, ny-1));
    iz = max(1, min(iz, nz-1));

    % 8 neighboring voxels
    c000 = img(ix  , iy  , iz  );
    c100 = img(ix+1, iy  , iz  );
    c010 = img(ix  , iy+1, iz  );
    c110 = img(ix+1, iy+1, iz  );
    c001 = img(ix  , iy  , iz+1);
    c101 = img(ix+1, iy  , iz+1);
    c011 = img(ix  , iy+1, iz+1);
    c111 = img(ix+1, iy+1, iz+1);

    % Trilinear interpolation formula
    val = ...
        c000*(1-fx)*(1-fy)*(1-fz) + ...
        c100*fx    *(1-fy)*(1-fz) + ...
        c010*(1-fx)*fy    *(1-fz) + ...
        c110*fx    *fy    *(1-fz) + ...
        c001*(1-fx)*(1-fy)*fz     + ...
        c101*fx    *(1-fy)*fz     + ...
        c011*(1-fx)*fy    *fz     + ...
        c111*fx    *fy    *fz;
end

%% ============================================================================
function [src, det0, du, dv] = computeDeltas(geo, i, angles)
% Compute full 3D geometry for CT projection
% Outputs:
%   S     - X-ray source position in normalized voxel space
%   det0  - Detector origin (corner pixel)
%   dU    - Direction vector per detector U pixel
%   dV    - Direction vector per detector V pixel

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

    % Initial detector corner point
    det0 = [ -(DSD - DSO);
             dU * (-nU/2 + 0.5);
             dV * ( nV/2 - 0.5) ];
    du = [0; dU; 0];
    dv = [0; 0; -dV];

    % Apply detector roll-pitch-yaw rotation
    det0 = rotateRPY(det0, dRoll, dPitch, dYaw);
    du   = rotateRPY(du,   dRoll, dPitch, dYaw);
    dv   = rotateRPY(dv,   dRoll, dPitch, dYaw);

    % Apply detector offset
    det0(2) = det0(2) + offDetU;
    det0(3) = det0(3) + offDetV;

    % Apply gantry ZYZ Euler rotation
    det0 = rotateZYZ(det0, angles(1,i), angles(2,i), angles(3,i));
    du   = rotateZYZ(du,   angles(1,i), angles(2,i), angles(3,i));
    dv   = rotateZYZ(dv,   angles(1,i), angles(2,i), angles(3,i));
    src  = rotateZYZ(src,  angles(1,i), angles(2,i), angles(3,i));

    % Apply volume origin offset
    det0(1) = det0(1) - offOrigX;
    det0(2) = det0(2) - offOrigY;
    det0(3) = det0(3) - offOrigZ;
    src(1)  = src(1)  - offOrigX;
    src(2)  = src(2)  - offOrigY;
    src(3)  = src(3)  - offOrigZ;

    % Shift coordinate system to align with voxel grid
    det0 = det0 + [geo.sVoxel(1,1)/2 - geo.dVoxel(1,1)/2;
                   geo.sVoxel(2,1)/2 - geo.dVoxel(2,1)/2;
                   geo.sVoxel(3,1)/2 - geo.dVoxel(3,1)/2];
    src = src + [geo.sVoxel(1,1)/2 - geo.dVoxel(1,1)/2;
                 geo.sVoxel(2,1)/2 - geo.dVoxel(2,1)/2;
                 geo.sVoxel(3,1)/2 - geo.dVoxel(3,1)/2];

    % Normalize to voxel units
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

%% ============================================================================
function p = rotateRPY(p, r, pch, y)
% Roll-Pitch-Yaw rotation matrix (detector rotation)

    R = [
        cos(r)*cos(pch),  cos(r)*sin(pch)*sin(y)-sin(r)*cos(y),  cos(r)*sin(pch)*cos(y)+sin(r)*sin(y);
        sin(r)*cos(pch),  sin(r)*sin(pch)*sin(y)+cos(r)*cos(y),  sin(r)*sin(pch)*cos(y)-cos(r)*sin(y);
        -sin(pch),        cos(pch)*sin(y),                       cos(pch)*cos(y)
    ];
    p = R * p;
end

%% ============================================================================
function p = rotateZYZ(p, a, t, ps)
% ZYZ Euler rotation matrix (gantry rotation)
    R = [
        cos(a)*cos(t)*cos(ps)-sin(a)*sin(ps),  -cos(a)*cos(t)*sin(ps)-sin(a)*cos(ps),  cos(a)*sin(t);
        sin(a)*cos(t)*cos(ps)+cos(a)*sin(ps),  -sin(a)*cos(t)*sin(ps)+cos(a)*cos(ps),  sin(a)*sin(t);
        -sin(t)*cos(ps),                        sin(t)*sin(ps),                         cos(t)
    ];
    p = R * p;
end

%% ============================================================================
function d = maxdistanceCuboid(geo, i)
% Compute safe starting distance to skip empty space before entering the volume

    maxCubX = (geo.nVoxel(1,1)/2 + abs(geo.offOrigin(1,i))/geo.dVoxel(1,1));
    maxCubY = (geo.nVoxel(2,1)/2 + abs(geo.offOrigin(2,i))/geo.dVoxel(2,1));
    maxCubZ = (geo.nVoxel(3,1)/2 + abs(geo.offOrigin(3,i))/geo.dVoxel(3,1));
    
    max_voxel_dim = max([geo.dVoxel(1,1), geo.dVoxel(2,1), geo.dVoxel(3,1)]);
    d = max(geo.DSO(i)/max_voxel_dim - sqrt(maxCubX^2 + maxCubY^2 + maxCubZ^2), 0.0);
end
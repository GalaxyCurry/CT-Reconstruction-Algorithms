function img_denoised = tvDenoise(img, lambda, maxIter)
% tvDenoise - 3D Total Variation (TV) denoising
%
% Inputs:
%   img        - 3D input image (single precision) [X, Y, Z]
%   lambda     - TV regularization parameter (denoising strength)
%   maxIter    - Maximum number of iterations
%
% Output:
%   img_denoised - Denoised 3D image

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% // update_u: 主变量更新
% u = (1-tau)*u + tau*( f + div(p)/lambda )
% 
% // update_p: 对偶变量更新
% q = p + tau2 * grad(u)
% p = q / max(1, ||q||)  // 投影L1 ball，保证TV约束
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% min u = λ⋅TV(u)+ 0.5∥u−f∥22
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------​
% Get image dimensions (X, Y, Z)
[X, Y, Z] = size(img);

% Voxel spacing 
dx = 1.0;
dy = 1.0;
dz = 1.0;

% Initialize variables
u = img;                            % Output denoised image
px = zeros(X, Y, Z, 'single');      % Dual variable for X-gradient
py = zeros(X, Y, Z, 'single');      % Dual variable for Y-gradient
pz = zeros(X, Y, Z, 'single');      % Dual variable for Z-gradient

% ====================== Main Iteration Loop ======================
for iter = 1:maxIter
    % Step-size parameters 
    tau2 = single(0.3 + 0.02 * (iter - 1));
    tau1 = single((1 / tau2) * (1/6 - 5 / (15 + (iter - 1))));
    
    % -------------------- Update u (primal variable) --------------------
    u_new = u;
    for z = 1:Z
        for y = 1:Y
            for x = 1:X
                % Compute divergence of the dual variables p
                div = computeDivergence(px, py, pz, x, y, z, dx, dy, dz);
                % Update rule for u
                u_new(x, y, z) = (1 - tau1) * u(x, y, z) + tau1 * (img(x, y, z) + div / lambda);
            end
        end
    end
    u = u_new;

    % -------------------- Update p (dual variables) --------------------
    for z = 1:Z
        for y = 1:Y
            for x = 1:X
                % Compute gradient of u
                [gx, gy, gz] = computeGradient(u, x, y, z, dx, dy, dz, X, Y, Z);
                % Compute unprojected dual variables
                qx = px(x, y, z) + tau2 * gx;
                qy = py(x, y, z) + tau2 * gy;
                qz = pz(x, y, z) + tau2 * gz;
                % Compute L2 norm for projection
                norm_q = sqrt(qx^2 + qy^2 + qz^2);
                norm_q = max(1.0, norm_q);
                % Project and update dual variables
                px(x, y, z) = qx / norm_q;
                py(x, y, z) = qy / norm_q;
                pz(x, y, z) = qz / norm_q;
            end
        end
    end
end

% Return final denoised image
img_denoised = u;
end

%% ===================== Divergence Computation =====================
function div = computeDivergence(px, py, pz, x, y, z, dx, dy, dz)
% Computes divergence of the vector field (px, py, pz) at voxel (x,y,z)
div = single(0.0);

% Z-component of divergence
if z > 1
    div = div + (pz(x, y, z) - pz(x, y, z-1)) / dz;
else
    div = div + pz(x, y, z);
end

% Y-component of divergence
if y > 1
    div = div + (py(x, y, z) - py(x, y-1, z)) / dy;
else
    div = div + py(x, y, z);
end

% X-component of divergence
if x > 1
    div = div + (px(x, y, z) - px(x-1, y, z)) / dx;
else
    div = div + px(x, y, z);
end
end

%% ===================== Gradient Computation =====================
function [gx, gy, gz] = computeGradient(u, x, y, z, dx, dy, dz, X, Y, Z)
% Computes 3D gradient of image u at voxel (x,y,z)
gx = single(0.0);
gy = single(0.0);
gz = single(0.0);

% X-direction gradient (forward difference)
if x < X
    gx = (u(x+1, y, z) - u(x, y, z)) / dx;
end

% Y-direction gradient (forward difference)
if y < Y
    gy = (u(x, y+1, z) - u(x, y, z)) / dy;
end

% Z-direction gradient (forward difference)
if z < Z
    gz = (u(x, y, z+1) - u(x, y, z)) / dz;
end
end
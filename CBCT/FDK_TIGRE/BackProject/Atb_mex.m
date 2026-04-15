function image = Atb_mex(projections, geo, alphas)
  
    
    nVoxelX = geo.nVoxel(1,1);
    nVoxelY = geo.nVoxel(2,1);
    nVoxelZ = geo.nVoxel(3,1);
    nDetecU = geo.nDetector(1,1);
    nDetecV = geo.nDetector(2,1);
    nalpha = length(alphas(1,:));
    
    image = zeros(nVoxelX, nVoxelY, nVoxelZ, 'single');
    projections = single(projections);

    for i = 1:nalpha
        fprintf('Processing projection %d/%d...\n', i, nalpha);
        
        alpha = alphas(1,i);
        theta = alphas(2,i);
        psi   = alphas(3,1);

        geo_curr.alpha = alpha;
        geo_curr.theta = theta;
        geo_curr.psi   = psi;

        geo_curr.dRoll = geo.rotDetector(1,i);
        geo_curr.dPitch = geo.rotDetector(2,i);
        geo_curr.dYaw = geo.rotDetector(3,i);

        geo_curr.offOrigX = geo.offOrigin(1,i);
        geo_curr.offOrigY = geo.offOrigin(2,i);
        geo_curr.offOrigZ = geo.offOrigin(3,i);

        geo_curr.offDetecU = geo.offDetector(1,i);
        geo_curr.offDetecV = geo.offDetector(2,i);

        geo_curr.DSD = geo.DSD(i);
        geo_curr.DSO = geo.DSO(i);
        geo_curr.COR = geo.COR(i);

        geo_curr.sVoxelX = geo.sVoxel(1,1);
        geo_curr.sVoxelY = geo.sVoxel(2,1);
        geo_curr.sVoxelZ = geo.sVoxel(3,1);

        geo_curr.dVoxelX = geo.dVoxel(1,1);
        geo_curr.dVoxelY = geo.dVoxel(2,1);
        geo_curr.dVoxelZ = geo.dVoxel(3,1);

        geo_curr.dDetecU = geo.dDetector(1,1);
        geo_curr.dDetecV = geo.dDetector(2,1);
        
      
        [xyzOrigin, deltaX, deltaY, deltaZ, S] = computeDeltasCube(geo_curr);
        
        % 预计算三角函数
        sinalpha = sin(alpha);
        cosalpha = cos(alpha);
        auxCOR = geo_curr.COR / geo_curr.dDetecU;
        
        proj = projections(:, :, i);
        
     
        for indZ = 1:nVoxelZ
            for indY = 1:nVoxelY
                for indX = 1:nVoxelX
                    
                    % ==============================
                    % 步骤1：计算旋转后的体素物理坐标
                    % ==============================
                    P.x = xyzOrigin.x + (indX-1)*deltaX.x + (indY-1)*deltaY.x + (indZ-1)*deltaZ.x;
                    P.y = xyzOrigin.y + (indX-1)*deltaX.y + (indY-1)*deltaY.y + (indZ-1)*deltaZ.y - auxCOR;
                    P.z = xyzOrigin.z + (indX-1)*deltaX.z + (indY-1)*deltaY.z + (indZ-1)*deltaZ.z;
                    
                    % ==============================
                    % 步骤2：计算射线向量（源 -> 体素）
                    % ==============================
                    vectX = P.x - S.x;
                    vectY = P.y - S.y;
                    vectZ = P.z - S.z;
                    
                    % ==============================
                    % 步骤3：射线与探测器平面求交
                    % ==============================
                    if abs(vectX) < 1e-12
                        continue; % 射线与X轴平行，无交点
                    end
                    t = (geo_curr.DSO - geo_curr.DSD - S.x) / vectX;
                    if t < 0
                        continue; % 交点在源后方，无效
                    end
                    y_proj = vectY * t + S.y;
                    z_proj = vectZ * t + S.z;
                    
                    % ==============================
                    % 步骤4：物理坐标转探测器像素坐标
                    % ==============================
                    u = y_proj + nDetecU * 0.5;
                    v = z_proj + nDetecV * 0.5;
                    
                    % 边界判断
                    if u < 1 || u >= nDetecU || v < 1 || v >= nDetecV
                        continue;
                    end
                    
                    % ==============================
                    % 步骤5：计算FDK反投影权重
                    % ==============================
                    % 计算真实体素坐标（未旋转）
                    realx = -(geo_curr.sVoxelX - geo_curr.dVoxelX)*0.5 + (indX-1)*geo_curr.dVoxelX + geo_curr.offOrigX;
                    realy = -(geo_curr.sVoxelY - geo_curr.dVoxelY)*0.5 + (indY-1)*geo_curr.dVoxelY + geo_curr.offOrigY + geo_curr.COR;
                    
                    
                    weight = (geo_curr.DSO + realy*sinalpha - realx*cosalpha) / geo_curr.DSO;
                    weight = 1 / (weight * weight);
                    
                    % ==============================
                    % 步骤6：双线性插值 + 反投影累加
                    % ==============================
                    val = bilinear_interp(proj, u, v);
                    image(indX, indY, indZ) = image(indX, indY, indZ) + val * weight;
                end
            end
        end
    end
end

%% 双线性插值
function val = bilinear_interp(img, u, v)
    [nV, nU] = size(img);
    
    % 插值坐标
    x0 = floor(u); y0 = floor(v);
    x1 = x0 + 1; y1 = y0 + 1;
    
    % 边界保护
    x0 = max(1, min(nU, x0)); x1 = max(1, min(nU, x1));
    y0 = max(1, min(nV, y0)); y1 = max(1, min(nV, y1));
    
    fx = u - x0; fy = v - y0;
    
    % 标准双线性插值
    val = (1-fx)*(1-fy)*img(y0, x0) ...
         + (1-fx)*fy*img(y0, x1) ...
         + fx*(1-fy)*img(y1, x0) ...
         + fx*fy*img(y1, x1);
end


%% ==============================
% computeDeltasCube
% 功能：计算旋转后的体素原点和增量
%% ==============================
function [xyzorigin, deltaX, deltaY, deltaZ, S] = computeDeltasCube(geo)
    % 初始化双精度点
    P.x = -(geo.sVoxelX/2 - geo.dVoxelX/2) + geo.offOrigX;
    P.y = -(geo.sVoxelY/2 - geo.dVoxelY/2) + geo.offOrigY;
    P.z = -(geo.sVoxelZ/2 - geo.dVoxelZ/2) + geo.offOrigZ;
    
    Px.x = P.x + geo.dVoxelX; Px.y = P.y;          Px.z = P.z;
    Py.x = P.x;          Py.y = P.y + geo.dVoxelY; Py.z = P.z;
    Pz.x = P.x;          Pz.y = P.y;          Pz.z = P.z + geo.dVoxelZ;
    
    % 物体ZYZT旋转
    P = eulerZYZT_matlab(geo, P);
    Px = eulerZYZT_matlab(geo, Px);
    Py = eulerZYZT_matlab(geo, Py);
    Pz = eulerZYZT_matlab(geo, Pz);
    
    % 探测器偏移
    P.z = P.z - geo.offDetecV;            P.y = P.y - geo.offDetecU;
    Px.z = Px.z - geo.offDetecV;          Px.y = Px.y - geo.offDetecU;
    Py.z = Py.z - geo.offDetecV;          Py.y = Py.y - geo.offDetecU;
    Pz.z = Pz.z - geo.offDetecV;          Pz.y = Pz.y - geo.offDetecU;
    
    % 探测器RPY旋转
    off = geo.DSO - geo.DSD;
    P.x = P.x + off;  Px.x = Px.x + off;  Py.x = Py.x + off;  Pz.x = Pz.x + off;
    P = rollPitchYawT_matlab(geo, P);
    Px = rollPitchYawT_matlab(geo, Px);
    Py = rollPitchYawT_matlab(geo, Py);
    Pz = rollPitchYawT_matlab(geo, Pz);
    P.x = P.x - off;  Px.x = Px.x - off;  Py.x = Py.x - off;  Pz.x = Pz.x - off;
    
    % 源点计算
    source.x = geo.DSD;
    source.y = -geo.offDetecU;
    source.z = -geo.offDetecV;
    source = rollPitchYawT_matlab(geo, source);
    source.x = source.x + off;
    
    % 缩放坐标到探测器像素单位
    P.z = P.z / geo.dDetecV;            P.y = P.y / geo.dDetecU;
    Px.z = Px.z / geo.dDetecV;          Px.y = Px.y / geo.dDetecU;
    Py.z = Py.z / geo.dDetecV;          Py.y = Py.y / geo.dDetecU;
    Pz.z = Pz.z / geo.dDetecV;          Pz.y = Pz.y / geo.dDetecU;
    source.z = source.z / geo.dDetecV;  source.y = source.y / geo.dDetecU;
    
    % 计算增量
    deltaX.x = Px.x - P.x; deltaX.y = Px.y - P.y; deltaX.z = Px.z - P.z;
    deltaY.x = Py.x - P.x; deltaY.y = Py.y - P.y; deltaY.z = Py.z - P.z;
    deltaZ.x = Pz.x - P.x; deltaZ.y = Pz.y - P.y; deltaZ.z = Pz.z - P.z;
    
    % 转单精度输出
    xyzorigin.x = single(P.x); xyzorigin.y = single(P.y); xyzorigin.z = single(P.z);
    S.x = single(source.x); S.y = single(source.y); S.z = single(source.z);
end

%% ==============================
% eulerZYZT旋转
%% ==============================
function p = eulerZYZT_matlab(geo, p)
    sin_alpha = sin(geo.alpha); cos_alpha = cos(geo.alpha);
    sin_theta = sin(geo.theta); cos_theta = cos(geo.theta);
    sin_psi = sin(geo.psi);   cos_psi = cos(geo.psi);
    
    x = p.x; y = p.y; z = p.z;
    
    p.x = x*(cos_psi*cos_theta*cos_alpha - sin_psi*sin_alpha) ...
        + y*(-cos_psi*cos_theta*sin_alpha - sin_psi*cos_alpha) ...
        + z*cos_psi*sin_theta;
    
    p.y = x*(sin_psi*cos_theta*cos_alpha + cos_psi*sin_alpha) ...
        + y*(-sin_psi*cos_theta*sin_alpha + cos_psi*cos_alpha) ...
        + z*sin_psi*sin_theta;
    
    p.z = -x*sin_theta*cos_alpha ...
        + y*sin_theta*sin_alpha ...
        + z*cos_theta;
end

%% ==============================
% rollPitchYawT旋转
%% ==============================
function p = rollPitchYawT_matlab(geo, p)
    sin_dRoll = sin(geo.dRoll);   cos_dRoll = cos(geo.dRoll);
    sin_dPitch = sin(geo.dPitch); cos_dPitch = cos(geo.dPitch);
    sin_dYaw = sin(geo.dYaw);     cos_dYaw = cos(geo.dYaw);
    
    x = p.x; y = p.y; z = p.z;
    
    p.x = cos_dRoll*cos_dPitch*x ...
        + sin_dRoll*cos_dPitch*y ...
        - sin_dPitch*z;
    
    p.y = (cos_dRoll*sin_dPitch*sin_dYaw - sin_dRoll*cos_dYaw)*x ...
        + (sin_dRoll*sin_dPitch*sin_dYaw + cos_dRoll*cos_dYaw)*y ...
        + cos_dPitch*sin_dYaw*z;
    
    p.z = (cos_dRoll*sin_dPitch*cos_dYaw + sin_dRoll*sin_dYaw)*x ...
        + (sin_dRoll*sin_dPitch*cos_dYaw - cos_dRoll*sin_dYaw)*y ...
        + cos_dPitch*cos_dYaw*z;
end
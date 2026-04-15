function volume = FDK(projections, angles, geo)

% 输入：
%   projections: 投影数据 [探测器行数, 探测器列数, 投影角度数]
%   angles: 投影角度数组 (弧度)
%   DSD: 源到探测器的距离 (mm)
%   DSO: 源到旋转中心的距离 (mm)
%   vol_size: 重建体积尺寸 [x, y, z]
%   voxel_size: 体素尺寸 [dx, dy, dz] (mm)
% 输出：
%   volume: 重建后的三维体积数据

DSD = geo.DSD(1,1);
DSO = geo.DSO(1,1);
vol_size = geo.nVoxel;
voxel_size = geo.dVoxel;


[num_rows, num_cols, num_angles] = size(projections);
vol_x = vol_size(1); vol_y = vol_size(2); vol_z = vol_size(3);
dx = voxel_size(1); dy = voxel_size(2); dz = voxel_size(3);

% 创建滤波器（斜坡滤波器加Hamming窗）
padded_length = 2*num_cols;
filter = design_filter(num_cols, padded_length);

% 初始化重建体积
volume = zeros(vol_x, vol_y, vol_z, 'single');

% 计算探测器参数
% 几何放大率：源到探测器 / 源到旋转中心
mag_factor = DSD / DSO;
% 探测器像素尺寸（匹配体素尺寸的放大）
du = voxel_size(1) * mag_factor; % 探测器x方向像素尺寸 (mm)
dv = voxel_size(2) * mag_factor; % 探测器y方向像素尺寸 (mm)
u_center = (num_cols + 1) / 2;
v_center = (num_rows + 1) / 2;

% 创建网格坐标
[x, y, z] = meshgrid(...
    ((0:vol_x-1) - (vol_x-1)/2) * dx, ...
    ((0:vol_y-1) - (vol_y-1)/2) * dy, ...
    ((0:vol_z-1) - (vol_z-1)/2) * dz);

% 主重建循环
for angle_idx = 1:num_angles
    beta = angles(angle_idx);
    fprintf('处理角度 %d/%d (%.1f°)\n', angle_idx, num_angles, rad2deg(beta));
    
    % 获取当前角度投影
    proj = projections(:, :, angle_idx);
    
    % 加权投影
    [U, V] = meshgrid(((1:num_cols) - u_center) * du, ((1:num_rows) - v_center) * dv);
    weight_factor = DSD ./ sqrt(DSD^2 + U.^2 + V.^2);
    weighted_proj = proj .* weight_factor;
    
    % 滤波处理
    filtered_proj = filter_projections(weighted_proj, filter, padded_length);
    
    % 三维反投影
    % 计算旋转矩阵
    R = [cos(beta)  sin(beta) 0; 
         -sin(beta) cos(beta) 0;
         0          0         1];
    
    % 坐标旋转
    coords = [x(:)'; y(:)'; z(:)'];
    rotated_coords = R * coords;
    xr = reshape(rotated_coords(1,:), size(x));
    yr = reshape(rotated_coords(2,:), size(y));
    zr = reshape(rotated_coords(3,:), size(z));
    
    % 计算投影坐标
    u_proj = (DSD * yr) ./ (DSO + xr) / du + u_center;
    v_proj = (DSD * zr) ./ (DSO + xr) / dv + v_center;
    
    % 反投影权重
    weight = DSD^2 ./ (DSO + xr).^2;
    
    % 插值获取投影值
    valid_idx = (u_proj >= 1) & (u_proj <= num_cols) & ...
                (v_proj >= 1) & (v_proj <= num_rows);
    
    proj_vals = zeros(size(x), 'single');
    proj_vals(valid_idx) = interp2(filtered_proj, u_proj(valid_idx), v_proj(valid_idx), 'linear', 0);
    
    % 累加到体积
    volume = volume + proj_vals .* weight;
end

% 角度积分归一化
volume = volume * (2 * pi / num_angles);
end


%%
function filter = design_filter(orig_length, padded_length)
% 频域斜坡滤波器（加Hamming窗）

% 生成频率轴
freq = linspace(-0.5, 0.5, orig_length);

% 斜坡滤波器（频域幅度与频率绝对值成正比）
ramlak = abs(freq);

% 加Hamming窗降低振铃效应
window = hamming(orig_length)';
filter_base = ramlak .* window;

% 初始化补零后的滤波器
filter = zeros(1, padded_length);  
half_orig = floor(orig_length/2);
half_padded = floor(padded_length/2);

% 将原始长度的滤波器居中放置到补零长度中
filter(half_padded - half_orig + 1 : half_padded + half_orig + mod(orig_length,2)) = filter_base;

% 转换为FFT后的频域顺序（移回直流分量到首位）
filter = fftshift(filter);

end

%%
function filtered_proj = filter_projections(proj, filter, padded_length)
% 对投影数据进行滤波
    [num_rows, num_cols] = size(proj);
    filtered_proj = zeros(num_rows, num_cols, 'single');
    
    % 扩展投影数据（避免边界效应）
    padded_proj = zeros(num_rows, padded_length);
    padded_proj(:, 1:num_cols) = proj;
    
    % 频域滤波：逐行处理
    for i = 1:num_rows
        row = padded_proj(i, :); 
        freq_row = fft(row);     
        filtered_freq = freq_row .* filter; 
        filtered_row = real(ifft(filtered_freq)); 
        % 截取原始长度
        filtered_proj(i, :) = filtered_row(1:num_cols); 
    end
end

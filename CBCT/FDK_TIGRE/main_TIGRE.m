% CBCT重建工程角度分析：多为"180+小扇角"有限角采集N张；
% CT分为有限角重建（多用于摆位受限）和稀疏角重建（多用于动态重建，心脏腹部等）

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Configuration & Switches
% =========================================================================
% --- Main control switch ---
% true: 投影数据代表"光子计数强度"(Photon Counting/Intensity)，需要做Log变换
% false: 投影数据代表"线积分"(Line Integral)，直接用于重建 
use_photon_counting_data = false; 

% use_file_input = true，从以下路径读取投影数据，否则运行模拟代码
use_file_input = false;
projection_file_path = 'path/to/your/projections.h5'; % (.h5, .mat, .raw等)

% --- Log transformation parameter ---
I0 = 12000; 
% =========================================================================


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Preparation

num_angles = 360; % 投影角度数
scan_angles = 2*pi; % 扫描角度
angles = linspace(0, scan_angles, num_angles);

DSD = 1160; 
DSO = 800; 
DSD = expandToLength(DSD, num_angles);
DSO = expandToLength(DSO, num_angles);

geo = defaultGeometry('nVoxel', [128,128,128], 'sVoxel', [32,32,32], 'nDetector', [128,128], 'sDetector', [64,64], 'mode', 'cone', 'DSD', DSD, 'DSO', DSO);
geo = checkGeo(geo, angles);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Module and Projection Simulation
phantom_vol = sheppLogan3D(geo.nVoxel,'Modified Shepp-Logan');

if use_file_input
   
    fprintf('Reading projection file: %s\n', projection_file_path);
    
    if exist(projection_file_path, 'file')
        load(projection_file_path, 'projections'); 
    else
        error('File does not exist: %s', projection_file_path);
    end
    
else
    fprintf('Simulating projection data...\n');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 输出的投影为积分图
%% Manual (Maybe have some wrong?????)
    % projections = simulate_projections(phantom_vol, angles, geo);
    
%% TIGRE "siddon" or "interpolated" or "JosephPlane"
    % projections = Ax(phantom_vol, geo, angles, 'siddon');
    % projections = Ax(phantom_vol, geo, angles, 'interpolated');
    projections = Ax(phantom_vol, geo, angles, 'JosephPlane');
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Log Transformation

if use_photon_counting_data
    fprintf('Photon counting intensity data detected, performing log transformation...\n');
    projections = intensityToLineIntegral(projections, I0);  
else
    fprintf('Data is in line integral format, skipping log transformation and proceeding directly to reconstruction.\n');
end







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Photon noise + electronic noise
if ~use_file_input % Add noise only in simulation mode to avoid repeatedly adding noise to real data
    noise_projections = addCTnoise(projections);
else
    noise_projections = projections;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDK Recon
% reconFDK = FDK_TIGRE(projections_josephplane,angles, geo, "filter", 'ram-lak');
reconFDKnoise = FDK_TIGRE(noise_projections, angles, geo, "filter", 'ram-lak');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Postprocessing
%% 3D deNoise
reconFDKnoise_denoised = im3DDenoise(reconFDKnoise,'TV',100,15);

%% cropCBCT: Crops all the sections that lie outide the area that is covered by the cone
% reconFDKnoise_denoised_croped = cropCBCT(reconFDKnoise_denoised,geo);









%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show 2D and 3D results
%% 2D results
figure;
subplot(1,5,1);
imshow(squeeze(phantom_vol(:,:,round(geo.nVoxel(3)/2))), []);
title('原始体模（中心切片）');
% subplot(1,5,2);
% imshow(squeeze(reconFDK(:,:,round(geo.nVoxel(3)/2))), []);
% title('无噪重建结果（中心切片）');
subplot(1,5,3);
imshow(squeeze(reconFDKnoise(:,:,round(geo.nVoxel(3)/2))), []);
title('有噪重建结果（中心切片）');
subplot(1,5,4);
imshow(squeeze(reconFDKnoise_denoised(:,:,round(geo.nVoxel(3)/2))), []);
title('有噪重建后去噪结果（中心切片）');
subplot(1,5,5);
imshow(squeeze(reconFDKnoise_denoised_croped(:,:,round(geo.nVoxel(3)/2))), []);
title('有噪重建后去噪且裁剪的结果（中心切片）');

%% 3D results
% plotImg(reconFDK,'Dim','Z');
















%% Length Extension Function
function [output] = expandToLength(input, targetLength)
    if length(input) == 1
        output = repmat(input, 1, targetLength);
    elseif length(input) == targetLength
        output = input;
    else
        error('Length mismatch! Must be 1 or %d', targetLength);
    end
end
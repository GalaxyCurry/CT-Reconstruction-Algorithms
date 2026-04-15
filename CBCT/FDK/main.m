% CBCT重建工程角度分析：多为"180+小扇角"有限角采集512张；
% CT分为有限角重建（多用于摆位受限）和稀疏角重建（多用于动态重建，心脏腹部等）

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Preparation
num_angles = 360; % 投影角度数
scan_angles = 2*pi; % 扫描角度
angles = linspace(0, scan_angles, num_angles);

DSD = 500; 
DSO = 200; 
DSD = expandToLength(DSD, num_angles);
DSO = expandToLength(DSO, num_angles);

geo = defaultGeometry('nVoxel', [128,128,128], 'sVoxel', [64,64,64], 'nDetector', [128,128], 'sDetector', [64,64], 'mode', 'cone', 'DSD', DSD, 'DSO', DSO);
geo = checkGeo(geo, angles);
geo.filter = 'ram-lak';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Module and Projection Simulation
% phantom_vol = sheppLogan3D(geo.nVoxel,'Modified Shepp-Logan');
% projections = simulate_projections(phantom_vol, angles, geo);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Photon noise + electronic noise
% noise_projections=addCTnoise(projections);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDK Recon
% reconFDK = FDK(projections, angles, geo);
% reconFDKnoise = FDK(noise_projections, angles, geo);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Postprocessing
%% 3D deNoise
% reconFDKnoise_denoised=im3DDenoise(reconFDKnoise,'TV',100,15);

%% cropCBCT: Crops all the sections that lie outide the area that is covered by the cone
reconFDKnoise_denoised_croped=cropCBCT(reconFDKnoise_denoised,geo);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show 2D and 3D results
figure;
subplot(1,5,1);
imshow(squeeze(phantom_vol(:,:,round(geo.nVoxel(3)/2))), []);
title('原始体模（中心切片）');
subplot(1,5,2);
imshow(squeeze(reconFDK(:,:,round(geo.nVoxel(3)/2))), []);
title('无噪重建结果（中心切片）');
subplot(1,5,3);
imshow(squeeze(reconFDKnoise(:,:,round(geo.nVoxel(3)/2))), []);
title('有噪重建结果（中心切片）');
subplot(1,5,4);
imshow(squeeze(reconFDKnoise_denoised(:,:,round(geo.nVoxel(3)/2))), []);
title('有噪重建后去噪结果（中心切片）');
subplot(1,5,5);
imshow(squeeze(reconFDKnoise_denoised_croped(:,:,round(geo.nVoxel(3)/2))), []);
title('有噪重建后去噪且裁剪的结果（中心切片）');

%plotImg(reconFDK,'Dim','Z');







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




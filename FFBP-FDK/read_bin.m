

filename = 'D:\tasks\CT_Recon\CBCT_FDK_GPU_Recon\ct-reconstruction-gpu-main\viz\vol_full.bin'; 
dataType = 'float';         % 根据实际数据类型修改 (常见: uint8/uint16/single)
byteOrder = 'ieee-le';       

%% 读取数据
imgss = readBinImage(filename, dataType, byteOrder);

%% 显示基本信息
fprintf('图像尺寸: %dx%dx%d\n', size(imgss,1), size(imgss,2), size(imgss,3));
fprintf('数据类型: %s\n', class(imgss));
fprintf('数值范围: [%.2f, %.2f]\n', min(imgss(:)), max(imgss(:)));


function img = readBinImage(filename, dataType, byteOrder)
    % readBinImage 读取 128x128x128 的二进制图像文件
    %
    % 输入参数:
    %   filename  - 二进制文件路径 (字符串)
    %   dataType  - 数据类型 (默认: 'uint16'，可选: 'uint8', 'int16', 'single', 'double' 等)
    %   byteOrder - 字节序 (默认: 'ieee-le' 小端序，可选: 'ieee-be' 大端序)
    %
    % 输出参数:
    %   img - 128x128x128 的三维图像数组
    
    % 设置默认参数
    if nargin < 2, dataType = 'uint16'; end
    if nargin < 3, byteOrder = 'ieee-le'; end
    
    % 打开文件 (rb 表示只读二进制模式)
    fID = fopen(filename, 'rb');
    if fID == -1
        error('无法打开文件，请检查路径: %s', filename);
    end
    
    % 定义图像尺寸
    imgSize = [128, 128, 128];
    numElements = prod(imgSize); % 计算总元素数: 128^3 = 2097152
    
    % 读取二进制数据
    % 语法: fread(文件ID, 元素数, 数据类型, 跳过字节, 字节序)
    img = fread(fID, numElements, dataType, 0, byteOrder);
    
    % 关闭文件
    fclose(fID);
    
    % 检查数据完整性
    if numel(img) ~= numElements
        warning('读取元素数 (%d) 与预期 (%d) 不符，请检查文件大小或数据类型', ...
            numel(img), numElements);
    end
    
    % 将一维数组重塑为 128x128x128 三维数组
    % MATLAB 按列优先 (Column-major) 存储，需与文件写入顺序一致
    img = reshape(img, imgSize);
end
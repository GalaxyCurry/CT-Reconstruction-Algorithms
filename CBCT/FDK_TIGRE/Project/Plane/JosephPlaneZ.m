function proj = JosephPlaneZ(img, geo, angles)

tic;

dx=geo.dVoxel(1);
dy=-geo.dVoxel(2);
dz=-geo.dVoxel(3);

Xplane=((0:geo.nVoxel(1))-size(img,1)/2)*dx; 
Yplane=((0:geo.nVoxel(2))-size(img,2)/2)*dy;
Zplane=((0:geo.nVoxel(3))-size(img,3)/2)*dz;
Xplane=Xplane-dx/2;
Yplane=Yplane-dy/2;
Zplane=Zplane-dz/2;


proj=zeros(geo.nDetector(1),geo.nDetector(2),size(angles,2));


for angle_index=1:size(angles,2)

    DSO=geo.DSO(1,angle_index);
    DSD=geo.DSD(1,angle_index);
    theta = angles(1, angle_index);  

    % ===================== 修改1：源位置（X轴正半轴，绕Z旋转）=====================
    % 初始位置(DSO, 0, 0)，绕Z轴旋转theta后的坐标
    SourceX = DSO * cos(theta);
    SourceY = DSO * sin(theta);
    SourceZ = 0;

    % ===================== 修改2：探测器中心位置（X轴负半轴，绕Z旋转）=====================
    % 初始位置(-(DSD-DSO), 0, 0)，绕Z轴旋转theta后的坐标
    DetectorX = -(DSD - DSO) * cos(theta);
    DetectorY = -(DSD - DSO) * sin(theta);
    DetectorZ = 0;

    % ===================== 修改3：探测器像素坐标计算（基于旋转方向向量）=====================
    DetectorLengthH=(floor(-geo.nDetector(1)/2):floor(geo.nDetector(1)/2))*geo.dDetector(1); 
    DetectorLengthV=(floor(-geo.nDetector(2)/2):floor(geo.nDetector(2)/2))*geo.dDetector(2); 

    % 探测器方向向量（绕Z轴旋转后保持正交）
    u_dir = [-sin(theta); cos(theta); 0];  % 水平方向（原Y轴方向旋转后）
    v_dir = [0; 0; 1];                      % 垂直方向（Z轴方向不变）

    % 生成探测器像素网格 [3, nu, nv]
    nu = length(DetectorLengthH);
    nv = length(DetectorLengthV);
    DetectorIndex = zeros(3, nu, nv);
    for iu = 1:nu
        for iv = 1:nv
            DetectorIndex(:, iu, iv) = [DetectorX; DetectorY; DetectorZ] + ...
                u_dir * DetectorLengthH(iu) + ...
                v_dir * DetectorLengthV(iv);
        end
    end

    % 保持与原代码一致的处理：裁剪最后一个元素（因DetectorLengthH长度为nu+1）
    DetectorIndex = DetectorIndex(:, 1:end-1, 1:end-1);
    NumberOfDetectorPixels(1) = size(DetectorIndex, 2);
    NumberOfDetectorPixels(2) = size(DetectorIndex, 3);

    % ===================== 修改4：探测器边界计算（基于方向向量）=====================
    % 水平方向（u）边界
    DetectorBoundaryU1 = DetectorIndex - u_dir * (geo.dDetector(1) / 2);
    DetectorBoundaryU2 = DetectorIndex + u_dir * (geo.dDetector(1) / 2);
    % 垂直方向（v）边界
    DetectorBoundaryV1 = DetectorIndex - v_dir * (geo.dDetector(2) / 2);
    DetectorBoundaryV2 = DetectorIndex + v_dir * (geo.dDetector(2) / 2);

    % ===================== 以下逻辑与原代码完全一致 =====================
    if(abs(SourceY-DetectorY)>=abs(SourceX-DetectorX) && abs(SourceY-DetectorY)>=abs(SourceZ-DetectorZ))
        SlopesU1=(SourceX-DetectorBoundaryU1(1,:,:))./(SourceY-DetectorBoundaryU1(2,:,:));
        InterceptU1=-SlopesU1*SourceY+SourceX;
        SlopesU2=(SourceX-DetectorBoundaryU2(1,:,:))./(SourceY-DetectorBoundaryU2(2,:,:));
        InterceptU2=-SlopesU2*SourceY+SourceX;
        SlopesV1=(SourceZ-DetectorBoundaryV1(3,:,:))./(SourceY-DetectorBoundaryV1(2,:,:));
        InterceptV1=-SlopesV1*SourceY+SourceZ;
        SlopesV2=(SourceZ-DetectorBoundaryV2(3,:,:))./(SourceY-DetectorBoundaryV2(2,:,:));
        InterceptV2=-SlopesV2*SourceY+SourceZ;
        intersection_slope1=(SourceX-DetectorIndex(1,:,:))./(SourceY-DetectorIndex(2,:,:));
        intersection_slope2=(SourceZ-DetectorIndex(3,:,:))./(SourceY-DetectorIndex(2,:,:));
        intersection_length=squeeze(abs(dy)./(cos(atan(intersection_slope1)).*cos(atan(intersection_slope2))));
        for iy=1:geo.nVoxel(2)
            coordX1=squeeze(SlopesU1)*(Yplane(iy)+dy/2)+squeeze(InterceptU1);
            coordX2=squeeze(SlopesU2)*(Yplane(iy)+dy/2)+squeeze(InterceptU2);
            coordZ1=squeeze(SlopesV1)*(Yplane(iy)+dy/2)+squeeze(InterceptV1);
            coordZ2=squeeze(SlopesV2)*(Yplane(iy)+dy/2)+squeeze(InterceptV2);
            image_x_index1=floor((coordX1-Xplane(1)+dx)./dx);
            image_x_index2=floor((coordX2-Xplane(1)+dx)./dx);
            image_z_index1=floor((coordZ1-Zplane(1)+dz)./dz);
            image_z_index2=floor((coordZ2-Zplane(1)+dz)./dz);
            proj(:,:,angle_index)=proj(:,:,angle_index)+Project_on_Y(img,...
                coordX1,coordX2,coordZ1,coordZ2,Xplane,Zplane,image_x_index1,...
                image_x_index2,image_z_index1,image_z_index2,dx,dz,iy).*...
                intersection_length;%./ray_normalization;
        end
    elseif(abs(SourceX-DetectorX)>=abs(SourceY-DetectorY) && abs(SourceX-DetectorX)>=abs(SourceZ-DetectorZ))
        SlopesU1=(SourceY-DetectorBoundaryU1(2,:,:))./(SourceX-DetectorBoundaryU1(1,:,:));
        InterceptU1=-SlopesU1*SourceX+SourceY;
        SlopesU2=(SourceY-DetectorBoundaryU2(2,:,:))./(SourceX-DetectorBoundaryU2(1,:,:));
        InterceptU2=-SlopesU2*SourceX+SourceY;
        SlopesV1=(SourceZ-DetectorBoundaryV1(3,:,:))./(SourceX-DetectorBoundaryV1(1,:,:));
        InterceptV1=-SlopesV1*SourceX+SourceZ;
        SlopesV2=(SourceZ-DetectorBoundaryV2(3,:,:))./(SourceX-DetectorBoundaryV2(1,:,:));
        InterceptV2=-SlopesV2*SourceX+SourceZ;
        intersection_slope1=(SourceY-DetectorIndex(2,:,:))./(SourceX-DetectorIndex(1,:,:));
        intersection_slope2=(SourceZ-DetectorIndex(3,:,:))./(SourceX-DetectorIndex(1,:,:));
        intersection_length=squeeze(abs(dx)./(cos(atan(intersection_slope1)).*cos(atan(intersection_slope2))));
        for ix=1:geo.nVoxel(1)
            coordY1=squeeze(SlopesU1)*(Xplane(ix)+dx/2)+squeeze(InterceptU1);
            coordY2=squeeze(SlopesU2)*(Xplane(ix)+dx/2)+squeeze(InterceptU2);
            coordZ1=squeeze(SlopesV1)*(Xplane(ix)+dx/2)+squeeze(InterceptV1);
            coordZ2=squeeze(SlopesV2)*(Xplane(ix)+dx/2)+squeeze(InterceptV2);
            image_y_index1=floor((coordY1-Yplane(1)+dy)./dy);
            image_y_index2=floor((coordY2-Yplane(1)+dy)./dy);
            image_z_index1=floor((coordZ1-Zplane(1)+dz)./dz);
            image_z_index2=floor((coordZ2-Zplane(1)+dz)./dz);
            proj(:,:,angle_index)=proj(:,:,angle_index)+Project_on_X(img,...
                coordY1,coordY2,coordZ1,coordZ2,Yplane,Zplane,image_y_index1,...
                image_y_index2,image_z_index1,image_z_index2,dx,dz,ix).*...
                intersection_length;%./ray_normalization;
        end
    else
        SlopesU1=(SourceX-DetectorBoundaryU1(1,:,:))./(SourceZ-DetectorBoundaryU1(3,:,:));
        InterceptU1=-SlopesU1*SourceZ+SourceX;
        SlopesU2=(SourceX-DetectorBoundaryU2(1,:,:))./(SourceZ-DetectorBoundaryU2(3,:,:));
        InterceptU2=-SlopesU2*SourceZ+SourceX;
        SlopesV1=(SourceY-DetectorBoundaryV1(2,:,:))./(SourceZ-DetectorBoundaryV1(3,:,:));
        InterceptV1=-SlopesV1*SourceZ+SourceY;
        SlopesV2=(SourceY-DetectorBoundaryV2(2,:,:))./(SourceZ-DetectorBoundaryV2(3,:,:));
        InterceptV2=-SlopesV2*SourceZ+SourceY;
        intersection_slope1=(SourceX-DetectorIndex(1,:,:))./(SourceZ-DetectorIndex(3,:,:));
        intersection_slope2=(SourceY-DetectorIndex(2,:,:))./(SourceZ-DetectorIndex(3,:,:));
        intersection_length=squeeze(abs(dz)./(cos(atan(intersection_slope1)).*cos(atan(intersection_slope2))));
        for iz=1:geo.nVoxel(3)
            coordX1=squeeze(SlopesU1)*(Zplane(iz)+dz/2)+squeeze(InterceptU1);
            coordX2=squeeze(SlopesU2)*(Zplane(iz)+dz/2)+squeeze(InterceptU2);
            coordY1=squeeze(SlopesV1)*(Zplane(iz)+dz/2)+squeeze(InterceptV1);
            coordY2=squeeze(SlopesV2)*(Zplane(iz)+dz/2)+squeeze(InterceptV2);
            image_x_index1=floor((coordX1-Xplane(1)+dx)./dx);
            image_x_index2=floor((coordX2-Xplane(1)+dx)./dx);
            image_y_index1=floor((coordY1-Yplane(1)+dy)./dy);
            image_y_index2=floor((coordY2-Yplane(1)+dy)./dy);
            proj(:,:,angle_index)=proj(:,:,angle_index)+Project_on_Z(img,...
                coordX1,coordX2,coordY1,coordY2,Xplane,Yplane,image_x_index1,...
                image_x_index2,image_y_index1,image_y_index2,dx,dy,iz).*...
                intersection_length;%./ray_normalization;
        end
    end
    fprintf('%d\n',angle_index);

end
toc
end

% ===================== 子函数100%原样保留，无任何修改 =====================
function proj=Project_on_Z(image,coordU1,coordU2,coordV1,coordV2,image_u_plane,...
    image_v_plane,image_index_u1,image_index_u2,image_index_v1,image_index_v2,...
    d_u_plane,d_v_plane,iz)
    proj=zeros(size(coordU1));
    tol_min=1e-6;
    for i=1:size(coordU1,1)
        for j=1:size(coordU1,2)
            p_value=0;
            s_index_u=min(image_index_u1(i,j),image_index_u2(i,j));
            e_index_u=max(image_index_u1(i,j),image_index_u2(i,j));
            s_index_v=min(image_index_v1(i,j),image_index_v2(i,j));
            e_index_v=max(image_index_v1(i,j),image_index_v2(i,j));
            for k=s_index_u:e_index_u
                if(k<1 || k>size(image,1))
                    continue;
                end
                if(s_index_u==e_index_u)
                    weight1=1;
                elseif(k==s_index_u)
                    weight1=(image_u_plane(k+1)-min(coordU1(i,j),coordU2(i,j)))/...
                        abs(coordU1(i,j)-coordU2(i,j));
                elseif(k==e_index_u)
                    weight1=(max(coordU1(i,j),coordU2(i,j))-image_u_plane(k))/...
                        abs(coordU1(i,j)-coordU2(i,j));
                else
                    weight1=abs(d_u_plane)/abs(coordU1(i,j)-coordU2(i,j));
                end
                if(abs(weight1)<tol_min)
                    weight1=0;
                end
                for l=s_index_v:e_index_v
                    if(l<1 || l>size(image,2))
                        continue;
                    end
                    if(s_index_v==e_index_v)
                        weight2=1;
                    elseif(l==s_index_v)
                        weight2=(max(coordV1(i,j),coordV2(i,j))-image_v_plane(l+1))/...
                            abs(coordV1(i,j)-coordV2(i,j));
                    elseif(l==e_index_v)
                        weight2=(image_v_plane(l)-min(coordV1(i,j),coordV2(i,j)))/...
                            abs(coordV1(i,j)-coordV2(i,j));
                    else
                        weight2=abs(d_v_plane)/abs(coordV1(i,j)-coordV2(i,j));
                    end
                    if(abs(weight2)<tol_min)
                        weight2=0;
                    end
                    p_value=p_value+(weight1*weight2)*image(l,k,iz);
                    assert( weight1>=0 && weight2>=0);
                end
            end
            proj(i,j)=p_value;
        end
    end
end

function proj=Project_on_X(image,coordU1,coordU2,coordV1,coordV2,image_u_plane,...
    image_v_plane,image_index_u1,image_index_u2,image_index_v1,image_index_v2,...
    d_u_plane,d_v_plane,ix)
    proj=zeros(size(coordU1));
    tol_min=1e-6;
    for i=1:size(coordU1,1)
        for j=1:size(coordU1,2)
            weight_sum=0;
            p_value=0;
            s_index_u=min(image_index_u1(i,j),image_index_u2(i,j));
            e_index_u=max(image_index_u1(i,j),image_index_u2(i,j));
            s_index_v=min(image_index_v1(i,j),image_index_v2(i,j));
            e_index_v=max(image_index_v1(i,j),image_index_v2(i,j));
            for k=s_index_u:e_index_u
                if(k<1 || k>size(image,2))
                    continue;
                end
                if(s_index_u==e_index_u)
                    weight1=1;
                elseif(k==s_index_u)
                    weight1=(max(coordU1(i,j),coordU2(i,j))-image_u_plane(k+1))/...
                        abs(coordU1(i,j)-coordU2(i,j));
                elseif(k==e_index_u)
                    weight1=(image_u_plane(k)-min(coordU1(i,j),coordU2(i,j)))/...
                        abs(coordU1(i,j)-coordU2(i,j));
                else
                    weight1=abs(d_u_plane)/abs(coordU1(i,j)-coordU2(i,j));
                end
                if(abs(weight1)<tol_min)
                    weight1=0;
                end
                for l=s_index_v:e_index_v
                    if(l<1 || l>size(image,3))
                        continue;
                    end
                    if(s_index_v==e_index_v)
                        weight2=1;
                    elseif(l==s_index_v)
                        weight2=(max(coordV1(i,j),coordV2(i,j))-image_v_plane(l+1))/...
                            abs(coordV1(i,j)-coordV2(i,j));
                    elseif(l==e_index_v)
                        weight2=(image_v_plane(l)-min(coordV1(i,j),coordV2(i,j)))/...
                            abs(coordV1(i,j)-coordV2(i,j));
                    else
                        weight2=abs(d_v_plane)/abs(coordV1(i,j)-coordV2(i,j));
                    end
                    if(abs(weight2)<tol_min)
                        weight2=0;
                    end
                    weight_sum=weight_sum+(weight1*weight2);
                    p_value=p_value+(weight1*weight2)*image(ix,k,l);
                end
            end
            proj(i,j)=p_value;
        end
    end
end
function proj=Project_on_Y(image,coordU1,coordU2,coordV1,coordV2,image_u_plane,...
    image_v_plane,image_index_u1,image_index_u2,image_index_v1,image_index_v2,...
    d_u_plane,d_v_plane,iy)

    proj=zeros(size(coordU1));
    tol_min=1e-6;
    for i=1:size(coordU1,1)
        for j=1:size(coordU1,2)
            p_value=0;
            s_index_u=min(image_index_u1(i,j),image_index_u2(i,j));
            e_index_u=max(image_index_u1(i,j),image_index_u2(i,j));
            s_index_v=min(image_index_v1(i,j),image_index_v2(i,j));
            e_index_v=max(image_index_v1(i,j),image_index_v2(i,j));
            for k=s_index_u:e_index_u
                if(k<1 || k>size(image,1))
                    continue;
                end
                if(s_index_u==e_index_u)
                    weight1=1;
                elseif(k==s_index_u)
                    weight1=(image_u_plane(k+1)-min(coordU1(i,j),coordU2(i,j)))/...
                        abs(coordU1(i,j)-coordU2(i,j));
                elseif(k==e_index_u)
                    weight1=(max(coordU1(i,j),coordU2(i,j))-image_u_plane(k))/...
                        abs(coordU1(i,j)-coordU2(i,j));
                else
                    weight1=abs(d_u_plane)/abs(coordU1(i,j)-coordU2(i,j));
                end
                if(abs(weight1)<tol_min)
                    weight1=0;
                end
                for l=s_index_v:e_index_v
                    if(l<1 || l>size(image,3))
                        continue;
                    end
                    if(s_index_v==e_index_v)
                        weight2=1;
                    elseif(l==s_index_v)
                        weight2=(max(coordV1(i,j),coordV2(i,j))-image_v_plane(l+1))/...
                            abs(coordV1(i,j)-coordV2(i,j));
                    elseif(l==e_index_v)
                        weight2=(image_v_plane(l)-min(coordV1(i,j),coordV2(i,j)))/...
                            abs(coordV1(i,j)-coordV2(i,j));
                    else
                        weight2=abs(d_v_plane)/abs(coordV1(i,j)-coordV2(i,j));
                    end
                    if(abs(weight2)<tol_min)
                        weight2=0;
                    end
                    p_value=p_value+(weight1*weight2)*image(k,iy,l);
                    assert(weight1>=0 && weight2>=0);
                    
                end
            end
            proj(i,j)=p_value;
        end
    end
end
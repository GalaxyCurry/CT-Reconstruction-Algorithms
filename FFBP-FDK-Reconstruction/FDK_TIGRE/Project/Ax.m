function [ projections ] = Ax(img, geo, angles, varargin )

% AX computes projections for images and geometry information
% 1.射线追踪
% 2.插值投影 
% 3.面积分加权


%% Optionals

ptype='Siddon';
% expected projection types. 'ray-voxel' is obsolete. Use 'Siddon'
expectedProjectionTypes = {'siddon','ray-voxel','interpolated','JosephPlane'};


if any(strcmp(varargin{1}, expectedProjectionTypes))
    ptype = varargin{1};
else
    assert(false,'Ax:InvalidInput','Projection type should be either ''interpolated'' or ''JosephPlane'' or ''Siddon''.');
end

%% image
assert(isa(img,'single'),'Ax:InvalidInput','Image should be single type');
assert(isreal(img),'Ax:InvalidInput','Image should be real (non-complex)');
assert(size(img,3)>1,'Ax:InvalidInput', 'image should be 3D');  
%% Angles
assert(isreal(angles),'Ax:InvalidInput','Angles should be real (non-complex)');
assert(size(angles,1)==1 | size(angles,1)==3 ,'Ax:InvalidInput','Angles should be of size 1xN or 3xN');
angles=double(angles); 
if size(angles,1)==1
   angles=repmat(angles,[3 1]);
   angles(2,:)=0;
   angles(3,:)=0;
end
%% geometry
assert(isequal([size(img,1) size(img,2) size(img,3)],squeeze(geo.nVoxel.')),'Ax:BadGeometry','nVoxel does not match with provided image size');

%% Temporary 

s= sum(abs(angles),2);
if (s(2)~=0 || s(3)~=0) && (strcmp(ptype,'Siddon')||strcmp(ptype,'ray-voxel')||strcmp(ptype,'JosephPlane')) && strcmp(geo.mode,'parallel')
    warning(['''Siddon'' Not supported for parallel beam arbitrary axis rotation, changed to ''JosephPlane'' or ''interpolated''.',newline,'Ignore this message if you are not purposedly triying enforce its use.']);
    ptype='interpolated';
end


%% TODO: 
switch lower(ptype)   
    case {'siddon', 'ray-voxel'}
        projections = Siddon(img, geo, angles);
    case 'interpolated'
        projections = Interpolated(img, geo, angles);
    case 'josephplane'
        projections = JosephPlaneZZ(img, geo, angles);
    otherwise
        error('Unknown projection type: %s. Supported types are: Siddon, ray-voxel, interpolated, JosephPlane', ptype);
end





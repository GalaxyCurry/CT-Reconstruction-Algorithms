function [ img ] = Atb(projections,geo,angles)
% ATB backprojection operator

%% image
assert(isreal(projections),'Atb:InvalidInput','Image should be real (non-complex)');
assert(size(projections,2)>1,'Atb:InvalidInput', 'Projections should be 2D');  
assert(size(projections,3)==size(angles,2),'Atb:InvalidInput', 'Number of projections should match number of angles.'); 
%% Angles
assert(isreal(angles),'Atb:InvalidInput','Angles should be real (non-complex)');
assert(size(angles,1)==1 | size(angles,1)==3 ,'Atb:InvalidInput','Angles should be of size 1xN or 3xN');
angles=double(angles); 
if size(angles,1)==1
   angles=repmat(angles,[3 1]);
   angles(2,:)=0;
   angles(3,:)=0;
end
%% geometry
assert(isequal([size(projections,2) size(projections,1)],geo.nDetector.'),'checkGeo:BadGeometry','nDetector does not match with provided image size');

%% backproject
img=Atb_mex(projections,geo,angles);

end


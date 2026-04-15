function [ imgO ] = im3DDenoise( img,type,varargin )
%   IMDENOISE3D removes noise of image with different methods
%   Currently only TV is supported. Input arguments are the image, the type
%   of denoising ('TV' only now) and the parameters for the denoising,
%   being number of iterations and hyperparameter currently available. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


if nargin>=4
    iter=varargin{1};
    hyper=varargin{2};
else
    iter=100;
    hyper=15;
end

if strcmp(type,'TV')
    % Normalize the image grayscale values to the [0, 1) range
    immin=min(img(:));
    img=(img-immin);
    immax=max(img(:));
    img=img./(immax+2*eps);

    % Generates an error if the data type of img is not single-precision float
    if ~isa(img, 'single')
        error('im3DDenoise: Input image of tvDenoise must be single precision');
    end
    if ismatrix(img)
        imgO=tvDenoise(cat(3,img,img),hyper,iter);
        imgO=imgO(:,:,1);
    else
        imgO=tvDenoise(img,hyper,iter);
    end
    clear img;
    
    imgO=imgO*immax;
    imgO=imgO+immin;
else
    error(['Unknown denoising method:  '  type ]);
end
clear img;

end




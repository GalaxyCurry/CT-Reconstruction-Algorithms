function projLine = intensityToLineIntegral(projection, I0, varargin)
% intensityToLineIntegral Converts CT detector photon intensity to line integral via logarithmic transformation
%
%   projLine = intensityToLineIntegral(projection)
%   projLine = intensityToLineIntegral(projection, I0)
%   projLine = intensityToLineIntegral(projection, I0, threshold)
%
% Description:
%   This function performs the standard pre-processing step in CT reconstruction.
%   It applies the Beer-Lambert law to convert measured X-ray intensity values
%   into the line integrals of the attenuation coefficient, which are required
%   by analytical reconstruction algorithms (FBP, FDK, etc.).
%
% Input Parameters:
%   projection - Input projection data (photon intensity I measured by the detector)
%                Can be a 2D matrix (single view) or 3D array (multiple views)
%   I0         - (Optional) Incident X-ray intensity (air scan / blank scan value)
%                Default: 12000
%   threshold  - (Optional) Lower bound for intensity values to avoid log(0)
%                Default: 1e-6
%
% Output Parameters:
%   projLine   - Output line integral data (i.e., -log(projection ./ I0))
%
% Physics:
%   Beer-Lambert Law: I = I0 * exp(-∫μ(x)dx)
%   Solved for line integral: ∫μ(x)dx = log(I0 / I)
%
% Example:
%   % Simulate intensity data
%   I0 = 12000;
%   projIntensity = I0 * exp(-rand(128, 128, 360));
%
%   % Convert to line integral using default parameters
%   projLine = intensityToLineIntegral(projIntensity);

% 1. Parse input arguments
if nargin < 2
    I0 = 12000; % Default incident intensity (Air scan value)
end

if nargin < 3
    threshold = 1e-6; % Default threshold for numerical stability
else
    threshold = varargin{1};
end

% 2. Apply threshold to prevent numerical issues (log(0) or division by zero)
val = max(projection, threshold);

% 3. Core calculation: Log transformation based on Beer-Lambert Law
projLine = log(I0 ./ val);

end
function [index_top, index_full] = calculate_tremor_index(data, faxis, varargin)
%CALCULATE_TREMOR_INDEX Calculate tremor index
%
% SYNTAX:
%   [index_top, index_full] = calculate_tremor_index(data, faxis)
%   [index_top, index_full] = calculate_tremor_index(data, faxis, ...
%                             'ArgumentName', Value, ...)
%
% INPUTS:
%   data    - Matrix with spectra (channel-by-frequency-by-time).
%   faxis   - Frequency axis in Hz (vector).
%
% OPTIONAL INPUT ARGUMENT-VALUE PAIRS:
%   'bandlim'   - Frequency band limits in Hz [f_low f_high]. Default is 
%                 [3.5 10].
%   'fraction'  - Fraction of index_full values to use for index_top.
%                 Default is 1/3.
%   'smoothing' - Smoothing paramaters passed on to the 'gausswin'
%                 function: gausswin(smoothing(1),smoothing(2)). Default is
%                 a window length of ~5 Hz and sigma = ~1 Hz.
%
% OUTPUTS:
%   index_top   - Scalar tremor index.
%   index_full  - Full tremor index vector.

% Check input:
narginchk(2,Inf)
validateattributes(data, {'single' 'double'}, {'3d'}, '', '''data''')
validateattributes(faxis, 'numeric', ...
    {'vector' 'numel' size(data,2) 'increasing'}, '', '''faxis''')

% Default values for optional input:
bandlim = [3.5 10];
fraction = 1/3;
smoothing = find(faxis-faxis(1)>5,1); % Window length ~5 Hz
smoothing(2) = smoothing/10;          % Sigma ~1 Hz

% Parse optional input:
if mod(length(varargin),2) % Check if the optional inputs come in pairs.
    error('Incomplete property-value pairs!');
else
    for i = 1:2:length(varargin) % Loop over pairs...
        switch lower(varargin{i})
            % Bandlimits
            case 'bandlim'
                bandlim = varargin{i+1};
                validateattributes(bandlim, 'numeric', {'vector' 'numel' 2 'increasing'}, ...
                    '', '''bandlim''')
            % Fraction
            case 'fraction'
                fraction = varargin{i+1};
                validateattributes(fraction, 'numeric', {'scalar' '>' 0 '<' 1}, ...
                    '', '''fraction''')
            % Smoothing
            case 'smoothing'
                smoothing = varargin{i+1};
                if ~isempty(smoothing)
                    validateattributes(smoothing, 'numeric', {'vector' 'numel' 2}, ...
                        '', '''smoothing''')
                end
        end
    end
end

% Dimension:
Nch = size(data,1);
Nomega = size(data,2);
Ntime = size(data,3);

% Smooth:
if ~isempty(smoothing)
    win = gausswin(smoothing(1),smoothing(2));
    win = win / sum(win(:));
    for iCh = 1:Nch
        data(iCh,:,:) = conv2(squeeze(data(iCh,:,:)), win, 'same');
    end
end

% Downsample:
ix = 2:2:Nomega;   % Skipping 0 Hz because it is usually not finite anyway...
faxis = faxis(ix);
Nomega = numel(faxis);
data = data(:,ix,:);
bFreq = faxis>bandlim(1) & faxis<bandlim(2);

% Index:
index_full = max(data(:,bFreq,:),[],2); % Max over frequency.
index_full = max(index_full,[],1); % Max over channel.
index_full = squeeze(index_full);
index_top = sort(index_full,1,'descend'); % Sort over time.
index_top = median(index_top(1:round(Ntime*fraction))); % Take the top 3rd.
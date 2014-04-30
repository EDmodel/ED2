function RGB = cmapping(varargin)
%CMAPPING   Colormap linear mapping/interpolation.
%
%   SYNTAX:
%           cmapping
%           cmapping(U)
%           cmapping(U,CMAP)
%           cmapping(U,CMAP,...,CNAN)
%           cmapping(U,CMAP,...,TYPE)
%           cmapping(U,CMAP,...,MODE)
%           cmapping(U,CMAP,...,MAPS)
%           cmapping(U,CMAP,...,CLIM)
%           cmapping(AX,...)
%     RGB = cmapping(...);
%
%   INPUT:
%     U     - May be one of the following options:
%              a) An scalar specifying the output M number of colors.
%              b) A vector of length M specifying the values at which
%                 the function CMAP(IMAP) will be mapped.
%              c) A matrix of size M-by-N specifying intensities to be
%                 mapped to an RGB (3-dim) image. May have NaNs elements. 
%             DEFAULT: Current colormap length.
%     CMAP  - A COLORMAP defined by its name or handle-function or RGB
%             matrix (with 3 columns) or by a combination of colors chars
%             specifiers ('kbcw', for example) to be mapped. See NOTE below
%             for more options.
%             DEFAULT: Current colormap
%     CNAN  - Color for NaNs values on U, specified by a 1-by-3 RGB color
%             or a char specifier.
%             DEFAULT: Current axes background (white color: [1 1 1])
%     TYPE  - String specifying the result type. One of:
%               'colormap'  Forces a RGB colormap matrix result (3 columns)
%               'image'     Forces a RGB image result (3 dimensions)
%             DEFAULT: 'image' if U is a matrix, otherwise is 'colormap'
%     MODE  - Defines the mapping way. One of:
%               'discrete'     For discrete colors
%               'continuous'   For continuous color (interpolates)
%             DEFAULT: 'continuous' (interpolates between colors)
%     MAPS  - Specifies the mapping type. One of (see NOTES below):
%               'scaled'   Scales mapping, also by using CLIM (as IMAGESC).
%               'direct'   Do not scales the mapping (as IMAGE).
%             DEFAULT: 'scaled' (uses CLIM)
%     CLIM  - Two element vector that, if given, scales the mapping within
%             this color limits. Ignored if 'direct' is specified.
%             DEFAULT: [0 size(CMAP,1)] or [0 1].
%     AX    - Uses specified axes or figure handle to set/get the colormap.
%             If used, must be the first input.
%             DEFAULT: gca
%
%   OUTPUT (all optional):
%     RGB - If U is not a matrix, this is an M-by-3 colormap matrix with
%           RGB colors in its rows, otherwise is an RGB image: M-by-N-by-3,
%           with the color red intensities defined by RGB(:,:,1), the green
%           ones by RGB(:,:,2) and the blue ones by RGB(:,:,3).
%
%   DESCRIPTION:
%     This functions has various functionalities like: colormap generator,
%     colormap expansion/contraction, color mapping/interpolation, matrix
%     intensities convertion to RGB image, etc.
%
%     The basic idea is a linear mapping between the CMAP columns
%     [red green blue] and the U data, ignoring its NaNs.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If a single value of U is required for interpolation, use [U U].
%     * If the char '-' is used before the CMAP name, the colors will be
%       flipped. The same occurs if U is a negative integer.
%
%   EXAMPLE:
%     % Colormaps:
%       figure, cmapping( 256,'krgby')            , colorbar
%       figure, cmapping(-256,'krgby' ,'discrete'), colorbar
%       figure, cmapping(log(1:100),[],'discrete'), colorbar
%     % Images:
%       u = random('chi2',2,20,30); u(15:16,7:9) = NaN;
%       u = peaks(30);  u(15:16,7:9) = NaN;
%       v = cmapping(u,jet(64),'discrete','k');
%       w = cmapping(u,cmapping(log(0:63),'jet','discrete'),'discrete');
%       figure, imagesc(u), cmapping(64,'jet'), colorbar
%        title('u')
%       figure, imagesc(v), cmapping(64,'jet'), colorbar
%        title('u transformed to RGB (look the colored NaNs)')
%       figure, imagesc(w) ,cmapping(64,'jet'), colorbar
%        title('u mapped with log(colormap)')
%       figure, imagesc(u), cmapping(log(0:63),'jet','discrete'), colorbar
%        title('u with log(colormap)')
%    
%   SEE ALSO:
%     COLORMAP, IND2RGB
%     and
%     CMJOIN by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cmapping.m
%   VERSION: 1.1 (Sep 02, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.0.1    Fixed little bug with 'm' magenta color. (Jun 30, 2009)
%   1.1      Fixed BUG with empty CMAP, thanks to Andrea Rumazza. (Sep 02,
%            2009) 

%   DISCLAIMER:
%   cmapping.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Sets defaults:
AX     = {};                     % Calculated inside.
U      = [];                     % Calculated inside.
CMAP   = [];                     % Calculated inside.
TYPE   = 'colormap';             % Changes to 'image' if U is a matrix.
CLIM   = [];                     % To use in scaling
CNAN   = [1 1 1];                % White 'w'
MODE   = 'continuous';           % Scaling to CLIM
MAPS   = 'scaled';               % Scaled mapping
method = 'linear';               % Interpolation method
mflip  = false;                  % Flip the colormap

% Gets figure handle and axes handle (just in case the default colormap or
% background color axes will be used.
HF     = get(0,'CurrentFigure');
HA     = [];
if ~isempty(HF)
 HA    = get(HF,'CurrentAxes');
 if ~isempty(HA)
  CNAN = get(HA,'Color');        % NaNs colors
 end
end

% Checks inputs:
if nargin>8
 error('CVARGAS:cmapping:tooManyInputs', ...
  'At most 8 inputs are allowed.')
elseif nargout>1
 error('CVARGAS:cmapping:tooManyOutputs', ...
  'At most 1 output is allowed.')
end

% Checks AX:
if (~isempty(varargin)) && ~isempty(varargin{1}) && ...
  (numel(varargin{1})==1) && ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'Type'),'axes')
 % Gets AX and moves all other inputs to the left:
 AX          = varargin(1);
 HA          = AX{1};
 CNAN        = get(HA,'Color');
 varargin(1) = [];
end

% Checks U:
Nargin = length(varargin);
if ~isempty(varargin)
 U           = varargin{1};
 varargin(1) = [];
end

% Checks CMAP:
if ~isempty(varargin)
 CMAP        = varargin{1};
 varargin(1) = []; 
end

% Checks input U, if not given uses as default colormap length:
% Note: it is not converted to a vector in case CMAP is a function and IMAP
%       was not given.
if isempty(U)
 % Gets default COLORMAP length:
 if ~isempty(HA)
  U = size(colormap(HA),1);
 else
  U = size(get(0,'DefaultFigureColormap'),1);
 end
elseif ndims(U)>2
 error('CVARGAS:cmapping:incorrectXInput', ...
  'U must be an scalar, a vector or a 2-dimensional matrix.')
end

% Checks input CMAP:
if isempty(CMAP)
 % CMAP empty, then uses default:
 if ~isempty(HA)
  CMAP = colormap(HA);
  if isempty(CMAP) % Fixed BUG, Sep 2009.
   CMAP = get(0,'DefaultFigureColormap');
   if isempty(CMAP)
    CMAP = jet(64);
   end
  end
 else
  CMAP = get(0,'DefaultFigureColormap');
  if isempty(CMAP)
   CMAP = jet(64);
  end
 end
 Ncmap = size(CMAP,1);
elseif isnumeric(CMAP)
 % CMAP as an [R G B] colormap:
 Ncmap = size(CMAP,1);
 if (size(CMAP,2)~=3) || ...
  ((min(CMAP(:))<0) || (max(CMAP(:))>1)) || any(~isfinite(CMAP(:)))
  error('CVARGAS:cmapping:incorrectCmapInput', ...
        'CMAP is an incorrect 3 columns RGB colors.')
 end
elseif ischar(CMAP)
 % String CMAP
 % Checks first character:
 switch CMAP(1)
  case '-'
   mflip = ~mflip;
   CMAP(1) = [];
   if isempty(CMAP)
    error('CVARGAS:cmapping:emptyCmapInput',...
     'CMAP function is empty.')
   end
 end
 if ~((exist(CMAP,'file')==2) || (exist(CMAP,'builtin')==5))
  % CMAP as a combination of color char specifiers:
  CMAP  = lower(CMAP);
  iy    = (CMAP=='y');
  im    = (CMAP=='m');
  ic    = (CMAP=='c');
  ir    = (CMAP=='r');
  ig    = (CMAP=='g');
  ib    = (CMAP=='b');
  iw    = (CMAP=='w');
  ik    = (CMAP=='k');
  Ncmap = length(CMAP);
  if (sum([iy im ic ir ig ib iw ik])~=Ncmap)
   error('CVARGAS:cmapping:incorrectCmapStringInput', ...
   ['String CMAP must be a valid colormap name or a combination of '...
    '''ymcrgbwk''.'])
  end
  % Convertion to [R G B]:
  CMAP       = zeros(Ncmap,3);
  CMAP(iy,:) = repmat([1 1 0],sum(iy),1);
  CMAP(im,:) = repmat([1 0 1],sum(im),1); % BUG fixed Jun 2009
  CMAP(ic,:) = repmat([0 1 1],sum(ic),1);
  CMAP(ir,:) = repmat([1 0 0],sum(ir),1);
  CMAP(ig,:) = repmat([0 1 0],sum(ig),1);
  CMAP(ib,:) = repmat([0 0 1],sum(ib),1);
  CMAP(iw,:) = repmat([1 1 1],sum(iw),1);
  CMAP(ik,:) = repmat([0 0 0],sum(ik),1);
 else
  % CMAP as a function name
  % Changes function name to handle:
  CMAP = str2func(CMAP);
  Ncmap = []; % This indicates a CMAP function input
 end
elseif isa(CMAP,'function_handle')
 Ncmap = []; % This indicates a CMAP function input
else
 % CMAP input unrecognized:
 error('CVARGAS:cmapping:incorrectCmapInput', ...
  'Not recognized CMAP input.') 
end

% Checks CMAP function handle:
if isempty(Ncmap)
 % Generates the COLORMAP from function:
 try
  temp = CMAP(2);
  if ~all(size(temp)==[2 3]) || any(~isfinite(temp(:))), error(''), end
  clear temp
 catch
  error('CVARGAS:cmapping:incorrectCmapFunction', ...
   ['CMAP function ''' func2str(CMAP) ''' must result in RGB colors.'])
 end
end

% Checks varargin:
while ~isempty(varargin)
 if     isempty(varargin{1})
  % continue
 elseif ischar(varargin{1})
  % string input
  switch lower(varargin{1})
   % CNAN:
   case 'y'         , CNAN = [1 1 0];
   case 'm'         , CNAN = [1 0 0];
   case 'c'         , CNAN = [0 1 1];
   case 'r'         , CNAN = [1 0 0];
   case 'g'         , CNAN = [0 1 0];
   case 'b'         , CNAN = [0 0 1];
   case 'w'         , CNAN = [1 1 1];
   case 'k'         , CNAN = [0 0 0];
   % MODE:
   case 'discrete'  , MODE = 'discrete';
   case 'continuous', MODE = 'continuous';
   % TYPE:
   case 'colormap'  , TYPE = 'colormap';
   case 'image'     , TYPE = 'image';
   % MAPS:
   case 'direct'    , MAPS = 'direct';
   case 'scaled'    , MAPS = 'scaled';
   % Incorrect input:
   otherwise
    error('CVARGAS:cmapping:incorrectStringInput',...
     ['Not recognized optional string input: ''' varargin{1} '''.'])
  end
 elseif isnumeric(varargin{1}) && all(isfinite(varargin{1}(:)))
  % numeric input
  nv = numel(varargin{1});
  if (nv==3) && (size(varargin{1},1)==1)
   % CNAN:
   CNAN = varargin{1}(:)';
   if (max(CNAN)>1) || (min(CNAN)<0)
    error('CVARGAS:cmapping:incorrectCnanInput',...
     'CNAN elements must be between 0 and 1.')
   end
  elseif (nv==2) && (size(varargin{1},1)==1)
   % CLIM:
   CLIM = sort(varargin{1},'ascend');
   if (diff(CLIM)==0)
    error('CVARGAS:cmapping:incorrectClimValues',...
     'CLIM must have 2 distint elements.')
   end
  else
   error('CVARGAS:cmapping:incorrectNumericInput',...
   'Not recognized numeric input.')
  end
 else
  error('CVARGAS:cmapping:incorrectInput',...
   'Not recognized input.')
 end
 % Clears current optional input:
 varargin(1) = [];
end % while


% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% U size:
[m,n] = size(U);
mn    = m*n;

% Checks TYPE:
if ~any([m n]==1)
 % Forces image TYPE if U is a matrix:
 TYPE = 'image';
elseif strcmp(TYPE,'colormap') && ~nargout && isempty(AX)
 % Changes the colormap on the specified or current axes if no output
 % argument:
 AX = {gca};
end

% Forces positive integer if U is an scalar, and flips CMAP if is negative:
if (mn==1)
 U = round(U);
 if (U==0)
  if ~nargout && strcmp(TYPE,'colormap')
   warning('CVARGAS:cmapping:incorrectUInput',...
    'U was zero and produces no colormap')
  else
   RGB = [];
  end
  return
 elseif (U<0)
  mflip = ~mflip;
  U     = abs(U);
 end
end

% Gets CMAP from function handle:
if isempty(Ncmap)
 if (mn==1)
  % From U:
  Ncmap = U(1);
 else
  % From default colormap:
  if ~isempty(HA)
   Ncmap = size(colormap(HA),1);
  else
   Ncmap = size(get(0,'DefaultFigureColormap'),1);
  end
 end
 CMAP = CMAP(Ncmap);
end

% Flips the colormap
if mflip
 CMAP = flipud(CMAP);
end

% Check CMAP when U is an scalar::
if (mn==1) && (U==Ncmap)
 % Finishes:
 if ~nargout && strcmp(TYPE,'colormap')
  if Nargin==0
   RGB = colormap(AX{:},CMAP);
  else
   colormap(AX{:},CMAP)
  end
 else
  RGB = CMAP;
  if strcmp(TYPE,'image')
   RGB = reshape(RGB,Ncmap,1,3);
  end
 end
 return
end

% Sets scales:
if strcmp(MAPS,'scaled')
 % Scaled mapping:
 if ~isempty(CLIM)
  if (mn==1)
   mn = U;
   U = linspace(CLIM(1),CLIM(2),mn)';
  else
   % Continue  
  end
 else
  CLIM = [0 1];
  if (mn==1)
   mn = U;
   U = linspace(CLIM(1),CLIM(2),mn)';
  else
   % Scales U to [0 1]:
   U = U-min(U(isfinite(U(:))));
   U = U/max(U(isfinite(U(:))));
   % Scales U to CLIM:
   U = U*diff(CLIM)+CLIM(1);
  end
 end
else
 % Direct mapping:
 CLIM = [1 Ncmap];
end

% Do not extrapolates:
U(U<CLIM(1)) = CLIM(1);
U(U>CLIM(2)) = CLIM(2);

% Sets CMAP argument:
umap = linspace(CLIM(1),CLIM(2),Ncmap)';

% Sets U:
if (mn==2) && (U(1)==U(2))
 % U = [Uo Uo] implicates U = Uo:
 U(2) = [];
 mn   = 1;
 m    = 1;
 n    = 1;
end

% Sets discretization:
if strcmp(MODE,'discrete')
 umap2 = linspace(umap(1),umap(end),Ncmap+1)';
 for k = 1:Ncmap
  U((U>umap2(k)) & (U<=umap2(k+1))) = umap(k);
 end
 clear umap2
end

% Forces column vector:
U = U(:);

% Gets finite data:
inan = ~isfinite(U);

% Initializes:
RGB  = repmat(reshape(CNAN,[1 1 3]),[mn 1 1]);

% Interpolates:
if (Ncmap>1) && (sum(~inan)>1)
 [Utemp,a,b]    = unique(U(~inan));
 RGBtemp = [...
  interp1(umap,CMAP(:,1),Utemp,method) ...
  interp1(umap,CMAP(:,2),Utemp,method) ...
  interp1(umap,CMAP(:,3),Utemp,method) ...
  ];
 RGB(~inan,:) = RGBtemp(b,:);
else
 % single color:
 RGB(~inan,1,:) = repmat(reshape(CMAP,[1 1 3]),[sum(~inan) 1 1]);
end

% Just in case
RGB(RGB>1) = 1; 
RGB(RGB<0) = 0;

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Output type:
if strcmp(TYPE,'colormap')
 RGB = reshape(RGB,mn,3);
 if ~isempty(AX)
  colormap(AX{:},RGB)
  if ~nargout 
   clear RGB
  end
 end
else
 RGB = reshape(RGB,[m n 3]);
end


% [EOF]   cmapping.m
function [CMAP,CLIM,WIDTH,REF,LEVELS] = ...
                                 cmfit(CMAP,CLIM,WIDTH,REF,CENTER,varargin) 
%CMFIT   Sets the COLORMAP and CAXIS to specific color bands. 
%
%   SYNTAX:
%                                      cmfit
%                                      cmfit(CMAP)
%                                      cmfit(CMAP,CLIM)
%                                      cmfit(CMAP,CLIM,WIDTH or LEVELS)
%                                      cmfit(CMAP,CLIM,WIDTH,REF)
%                                      cmfit(CMAP,CLIM,WIDTH,REF,CENTER)
%                                      cmfit(AX,...)
%     [CMAPF,CLIMF,WIDTHF,REFF,LEVF] = cmfit(...);
%
%   INPUT:
%     CMAP   - Fits the specified colormap function or RGB colors. 
%              DEFAULT: (current figure colormap)
%     CLIM   - 2 element vector spacifying the limits of CMAP. 
%              DEFAULT: (limits of a COLORBAR)
%     WIDTH  - Color band width (limits are computed with CAXIS) for each
%     or       band or a row vector specifying the LEVELS on each band (see
%     LEVELS   NOTE below).
%              DEFAULT: (fills the ticks of a COLORBAR)
%     REF    - Reference level to start any of the color bands.
%              DEFAULT: (generally the middle of CLIM)
%     CENTER - Logical specifying weather the colormap should be center in
%              the REF value or not.
%              DEFAULT: false (do not centers)
%     AX     - Uses the specified figure or axes handle.
%
%   OUTPUT (all optional):
%     CMAPF  - RGB fitted color map (with 3 columns).
%     CLIMF  - Limits of CMAPF.
%     WIDTHF - Width of fitted colorbands.
%     REFF   - Reference of fitted colorbands.
%     LEVF   - Levels for the color bands.
%
%   DESCRIPTION:
%     This program sets the current figure colormap with specified
%     band-widths of colors taking the CAXIS limits as reference. When the 
%     optional input argument CENTER is true, the colormap is moved and
%     expanded so its middle color will be right at REF. This will help for
%     distinguish between positive and negative values (REF=0).
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * When one of the first two inputs is missing, they are automatically
%       calculated by using a COLORBAR (created temporarly if necesary). In
%       this case CBHANDLE is necesary.
%     * When CMAP is used as output, the current figure colormap won't be
%       modificated. Use 
%         >> colormap(CMAP)
%       after this function, if necesary.
%     * When LEVELS are used instead of band WINDTH, it shoud be
%       monotonically increasing free of NaNs and of length equal to the
%       number of colors minus one, on the output colormap.
% 
%   SEE ALSO:
%     COLORMAP
%     and 
%     CMAPPING, CBFIT by Caros Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cmfit.m
%   VERSION: 1.0 (Jun 08, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)

%   DISCLAIMER:
%   cmfit.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2008,2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Sets defaults: 
AX  = {};    % Axes input
tol = 1;     % Adds this tolerance to the decimal precision
hfig = {get(0,'CurrentFigure')};
if ~isempty(hfig)
 hax = {get(hfig{1},'CurrentAxes')};
 if isempty(hax{1}), hax = {}; end
else
 hfig = {};
 hax  = {};
end

% Checks inputs:
if nargin>6
 error('CVARGAS:cmfit:tooManyInputs', ...
  'At most 6 inputs are allowed.')
end
if nargin>5
 error('CVARGAS:cmfit:tooManyOutputs', ...
  'At most 5 outputs are allowed.')
end

% Saves number of arguments:
Nargin = nargin;

% Checks AX input:
if (Nargin>0) && ~isempty(CMAP) && (numel(CMAP)==1) && ...
  ishandle(CMAP)
 % Gets AX and moves all other inputs to the left:
 AX = {CMAP};
 switch get(AX{1},'Type')
  case 'axes'
   hax  = AX;
   hfig = {get(hax{1},'Parent')};
  case {'figure','uipanel'}
   hfig = {AX{1}};
   hax  = {get(hfig{1},'CurrentAxes')};
   if isempty(hax{1}), hax = {}; end
  otherwise
   error('CVARGAS:cmfit:incorrectAxHandle',...
    'AX must be a valid axes or figure handle.')
 end
 if (Nargin>1)
  CMAP = CLIM;
  if (Nargin>2)
   CLIM = WIDTH;
   if (Nargin>3)
    WIDTH = REF;
    if (Nargin>4)
     REF = CENTER;
     if (Nargin>5)
      CENTER = varargin{1};
     end
    end
   end
  end
 end
 Nargin = Nargin-1;
end

% Checks CMAP input:
if Nargin<1 || isempty(CMAP)
 if ~isempty(hax)
  CMAP = colormap(hax{1});
 else
  CMAP = get(0,'DefaultFigureColormap');
 end
end

% Checks CLIM input:
if Nargin<2
 CLIM = [];
end

% Checks WIDTH input:
if Nargin<3
 WIDTH = [];
end

% Checks REf input:
if Nargin<4
 REF = [];
end

% Checks CENTER input:
if Nargin<5 || isempty(CENTER)
 CENTER = false;
end

% Look for WIDTH and REF from a (temporarly) colorbar:
if isempty(WIDTH) || (length(WIDTH)==1 && (isempty(REF) || ...
  (isempty(CLIM) && (isempty(hax) || ...
  ~strcmp(get(hax{1},'CLimMode'),'manual')))))
 if ~isempty(CLIM)
  caxis(hax{:},CLIM)
 end
 if ~isempty(AX) && ~isempty(cbhandle(AX{1}))
  h = cbhandle(AX{1}); doclear = false; h = h(1);
 elseif ~isempty(hax) && ~isempty(cbhandle(hax{1}))
  h = cbhandle(hax{1}); doclear = false; h = h(1);
 elseif ~isempty(hfig) && ~isempty(cbhandle(hfig{1}))
  h = cbhandle(hfig{1}); doclear = false; h = h(1);
 else
  h = colorbar; doclear = true;
 end
 ticks = get(h,'XTick');
 lim   = get(h,'XLim');
 if isempty(ticks)
  ticks = get(h,'YTick');
  lim   = get(h,'YLim');
 end
 if isempty(WIDTH)
  WIDTH = diff(ticks(1:2));
 end
 if isempty(CLIM)
  CLIM = lim;
 end
 if isempty(REF) && ~CENTER
  REF = ticks(1);
 end
 if doclear
  delete(h)
 end
end

% Centers at the middle:
if CENTER && isempty(REF)
 REF = 0;
end 

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------


% Gets minimum width from specified levels:
NL = length(WIDTH); 
if (NL>1)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NONLINEAR CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 LEVELS = sort(WIDTH);
 
 % Gets LEVELS width:
 wLEVELS = diff(LEVELS);
 
 % Scales to CLIM:
 if ~isempty(CLIM)
  % Scales to [0 1]
  LEVELS = LEVELS-LEVELS(1);
  LEVELS = LEVELS/LEVELS(end);
  % Scales to CLIM:
  LEVELS = LEVELS*diff(CLIM)+CLIM(1);
 else
  CLIM  = [LEVELS(1) LEVELS(end)]; 
 end
 
 % Gets precision:
 if isinteger(wLEVELS) % Allows integer input: uint8, etc. 
  wLEVELS = double(wLEVELS);
 end
 temp = warning('off','MATLAB:log:logOfZero');
 precision = floor(log10(abs(wLEVELS))); % wLEVELS = Str.XXX x 10^precision.
 precision(wLEVELS==0) = 0; % M=0 if x=0.
 warning(temp.state,'MATLAB:log:logOfZero')
 precision = min(precision)-tol;
 
 % Sets levels up to precision:
 wLEVELS = round(wLEVELS*10^(-precision));
 
 % Gets COLORMAP for each LEVEL:
 if CENTER
  % Centers the colormap:
  ind = (REF==LEVELS);
  if ~any(ind)
   error('CVARGAS:cmfit:uncorrectRefLevel',...
    'When CENTER, REF level must be on of the specifyied LEVELS.')
  end
  Nl     = sum(~ind(1:find(ind)));
  [Nl,l] = max([Nl (NL-1-Nl)]);
  wCMAP  = cmapping(2*Nl,CMAP);
  if l==1
   wCMAP = wCMAP(1:NL-1,:);
  else
   wCMAP = wCMAP(end-NL+2:end,:);
  end
 else
  wCMAP  = cmapping(NL-1,CMAP);
 end
 
 % Gets minimum band width:
 WIDTH = wLEVELS(1);
 for k = 1:NL-1
  wlev    = wLEVELS;
  wlev(k) = [];
  WIDTH   = min(min(gcd(wLEVELS(k),wlev)),WIDTH);
 end
 
 % Gets number of bands:
 wLEVELS = wLEVELS/WIDTH;
 
 % Gets new CMAP:
 N = sum(wLEVELS);
 try
  CMAP = repmat(wCMAP(1,:),N,1);
 catch
  error('CVARGAS:cmfit:memoryError',...
   ['The number of colors (N=' int2str(N) ') for the new colormap ' ...
    'is extremely large. Try other LEVELS.'])
 end
 ko = wLEVELS(1);
 for k = 2:NL-1;
  CMAP(ko+(1:wLEVELS(k)),:) = repmat(wCMAP(k,:),wLEVELS(k),1);
  ko = ko+wLEVELS(k);
 end
 
else
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % LINEAR CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Gets CLIM:
 if isempty(CLIM)
  CLIM = caxis(hax{:});
 end

 % Sets color limits to be a multipler of WIDTH passing through REF:
 N1   = ceil((+REF-CLIM(1))/WIDTH);
 N2   = ceil((-REF+CLIM(2))/WIDTH);
 CLIM = REF + [-N1 N2]*WIDTH;

 % Sets colormap with NC bands:
 Nc = round(diff(CLIM)/WIDTH);
 if CENTER
  % Necesary colorbands number to be centered:
  Nmin        = [N1 N2];
  [Nmax,imax] = max(Nmin);
  Nmin(imax)  = [];
  Nc2         = Nc + Nmax - Nmin;
  % Generate a colormap with this size:
  CMAP = cmapping(Nc2,CMAP);
  if imax==1
   CMAP = CMAP(1:Nc,:);
  else
   CMAP = flipud(CMAP);
   CMAP = CMAP(1:Nc,:);
   CMAP = flipud(CMAP);
  end
 else
  CMAP = cmapping(Nc,CMAP);
 end
 
 % Sets levels:
 LEVELS = linspace(CLIM(1),CLIM(2),size(CMAP,1))';
end

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------
if ~isempty(AX)
 colormap(AX{:},CMAP)
 caxis(AX{:},CLIM)
end
if ~nargout
 if isempty(AX)
  colormap(CMAP)
  caxis(CLIM)
 end
 clear CMAP
end


% [EOF]   cmfit.m
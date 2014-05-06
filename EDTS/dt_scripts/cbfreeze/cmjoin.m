function [CMAP,LEV,WID,CLIM] = cmjoin(varargin)
%CMJOIN   Joins colormaps at certain levels.
%
%   SYNTAX:
%                           cmjoin(CMAPS)
%                           cmjoin(CMAPS,LEV)
%                           cmjoin(CMAPS,LEV,WID)
%                           cmjoin(CMAPS,LEV,WID,CLIM)
%                           cmjoin(AX,...)
%     [CMAP,LEV,WID,CLIM] = cmjoin(...);
%
%   INPUT:
%     CMAPS - Cell with the N colormaps handles, names or RGB colors to be
%             joined. See NOTE below.
%     LEV   - One of:
%               a) N-1 scalars specifying the color levels where the
%                  colormaps will be joined (uses CAXIS). See NOTE below.
%               b) N integers specifying the number of colors for each
%                  colormap.
%               c) N+1 scalars specifying the color limits for each
%                  colormap (sets CAXIS). See NOTE below.
%             DEFAULT: Tries to generate a CMAP with default length.
%     WID   - May be one (or N) positive scalar specifying the width for
%             every (or each) color band. See NOTE below.
%             DEFAULT: uses CAXIS and LEV to estimate it.  
%     CLIM  - 2 elements row vector specifying the color limits values. May
%             be changed at the end, because of the discretization of the
%             colormaps.
%             DEFAULT: uses CAXIS or [0 1] if there are no axes.
%     AX    - Uses the specified axes handle to get/set the CMAPS. If used,
%             must be the first input.
%             DEFAULT: gca
%
%   OUTPUT (all optional):
%     CMAP - RGB colormap output matrix, M-by-3.
%     LEV  - Final levels used.
%     WID  - Final widths used.
%     CLIM - Final color limits used.
%
%   DESCRIPTION:
%     This function join two colormaps at specific level. Useful for
%     joining colormaps at zero (for example) and distinguish from positive
%     and negative values.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If no output is required or an axes handle were given, the current
%       COLORMAP and CAXIS are changed.
%     * If any of the inputs on CMAPS is a function name, 'jet', for
%       example, it can be used backwards (because CMAPPING is used) if
%       added a '-' at the beggining of its name: '-jet'.
%     * When LEV is type b) and WID is specifyed, the latter is taken as
%       relative colorbans widths between colormaps.
%
%   EXAMPLE:
%     figure(1), clf, surf(peaks)
%     cmjoin({'copper','-summer'},2.5)
%      shading interp, colorbar, axis tight, zlabel('Meters')
%      title('Union at 2.5 m')
%     %
%     figure(2), clf, surf(peaks) 
%     cmjoin({'copper','-summer'},2.5,0.5)
%      shading interp, colorbar, axis tight, zlabel('Meters')
%      title('Union at 2.5 m and different color for each 0.5 m band')
%     %
%     figure(3), clf, surf(peaks)
%     cmjoin({'copper','summer'},2.5,[2 0.5])
%      shading interp, colorbar, axis tight, zlabel('Metros')
%      title('Union at 2.5 m with lengths 2 and 0.5')
%     %
%     figure(4), clf, surf(peaks)
%     cmjoin({'copper','summer'},[-10 2.5 10],[2 0.5])
%      shading interp, colorbar, axis tight, zlabel('Metros')
%      title('Union at 2.5 m with lengths 2 and 0.5 and specified levels')
%     %
%     figure(5), clf, surf(peaks)
%     cmjoin({'copper','summer'},[10 8],[4 1])
%      shading interp, colorbar, axis tight, zlabel('Metros')
%      title('Union at 2.5 m with specified levels number of colors and widths 4:1')
%    
%   SEE ALSO:
%     COLORMAP, COLORMAPEDITOR
%     and
%     CMAPPING by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cmjoin.m
%   VERSION: 2.0 (Jun 08, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0       Released as SETCOLORMAP. (Nov 07, 2006)
%   1.1       English translation. (Nov 11, 2006)
%   2.0       Rewritten and renamed code (from SETCOLORMAPS to CMJOIN. Now
%             joins multiple colormaps. Inputs changed. (Jun 08, 2009)

%   DISCLAIMER:
%   cmjoin.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2006,2009 Carlos Adrian Vargas Aguilera

% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
tol  = 1;            % When rounding the levels.

% Checks inputs and outputs number:
if nargin<1
 error('CVARGAS:cmjoin:notEnoughInputs',...
  'At least 1 input is required.')
end
if nargin>5
 error('CVARGAS:cmjoin:tooManyInputs',...
  'At most 5 inputs are allowed.')
end
if nargout>4
 error('CVARGAS:cmjoin:tooManyOutputs',...
  'At most 4 outputs are allowed.')
end

% Checks AX:
AX = {get(get(0,'CurrentFigure'),'CurrentAxes')};
if isempty(AX{1})
 AX = {};
end
if (length(varargin{1})==1) && ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'Type'),'axes')
 AX = varargin(1);
 varargin(1) = [];
 if isempty(varargin)
  error('CVARGAS:cmjoin:notEnoughInputs',...
   'CMAPS input must be given.')
 end
end

% Checks CMAPS:
CMAPS  = varargin{1};
Ncmaps = length(CMAPS);
if ~iscell(CMAPS) || (Ncmaps<2)
 error('CVARGAS:cmjoin:incorrectCmapsType',...
  'CMAPS must be a cell input with at least 2 colormaps.')
end
varargin(1) = [];
Nopt        = length(varargin);

% Checks LEV and sets Ncol and Jlev:
Ncol = []; % Number of colors for each colormap.
Jlev = []; % Join levels.
LEV  = []; % Levels at which each CMAPS begins and ends.
if (Nopt<1) || isempty(varargin{1})
 % continue as empty
elseif ~all(isfinite(varargin{1}(:)))
 error('CVARGAS:cmjoin:incorrectLevValue',...
  'LEV must be integers or scalars.')
else
 Nopt1 = length(varargin{1}(:));
 if (Nopt1==Ncmaps)
  % Specifies number of colors:
  Ncol = varargin{1}(:);
  if ~all(Ncol==round(Ncol))
   error('CVARGAS:cmjoin:incorrectLevInput',...
    'LEV must be integers when defines number of colors.')
  end
 elseif ~all(sort(varargin{1})==varargin{1})
  error('CVARGAS:cmjoin:incorrectLevInput',...
   'LEV must be monotonically increasing.')
 elseif Nopt1==(Ncmaps-1) 
  Jlev = varargin{1}(:);
 elseif Nopt1==(Ncmaps+1)
  LEV = varargin{1}(:);
 else
  error('CVARGAS:cmjoin:incorrectLevLength',...
   'LEV must have any of length(CMAPS)+[-1 0 1] elements.')
 end
end

% Checks WID:
Tcol = []; % Total number of colors for output colormap.
if (Nopt<2) || isempty(varargin{2})
 % Tries to generate a colormap with default length with every colorband
 % of the same width:
 WID = [];
 if ~isempty(AX)
  Tcol = size(colormap(AX{:}),1);
 else
  Tcol = size(get(0,'DefaultFigureColormap'),1);
 end
else
 WID = varargin{2}(:);
 WID(~isfinite(WID) | (WID<0)) = 0;
 if ~any(WID>0)
  error('CVARGAS:cmjoin:incorrectWidInput',...
   'At least one WID must be positive.')
 end
 if length(WID)==1
  WID = repmat(abs(varargin{2}),Ncmaps,1);
 elseif length(WID)~=Ncmaps
  error('CVARGAS:cmjoin:incorrectWidLength',...
   'WID must have length 1 or same as CMAPS.')
 end
end

% Checks CLIM:
if (Nopt<3) || isempty(varargin{3})
 % Sets default CLIM:
 if ~isempty(LEV)
  CLIM = [LEV(1) LEV(end)];
 elseif ~isempty(AX)
  CLIM = caxis(AX{:});
 else
  CLIM = [0 1];
 end
else
 CLIM = varargin{3}(:).';
 if (length(CLIM)==2) && (diff(CLIM)>0) && isfinite(diff(CLIM))
  % continue
 else
  error('CVARGAS:cmjoin:incorrectClimInput',...
   'CLIM must be a valid color limits. See CAXIS for details.')
 end
end


% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Gets rounding precision:
temp = warning('off','MATLAB:log:logOfZero');
if ~isempty(WID)
 tempp               = WID;
 precision           = floor(log10(abs(tempp)));
 precision(tempp==0) = 0;
 precision           = min(precision)-tol;
 % Rounds:
 WID   = round(WID*10^(-precision))*10^precision;
 if ~isempty(LEV)
  LEV(1)        = floor(LEV(1)       *10^(-precision))*10^precision;
  LEV(2:end-1)  = round(LEV(2:end-1) *10^(-precision))*10^precision;
  LEV(end)      = ceil(LEV(end)      *10^(-precision))*10^precision;
 elseif ~isempty(Jlev)
  Jlev(1)       = floor(Jlev(1)      *10^(-precision))*10^precision;
  Jlev(2:end-1) = round(Jlev(2:end-1)*10^(-precision))*10^precision;
  Jlev(end)     = ceil(Jlev(end)     *10^(-precision))*10^precision;
 end
elseif ~isempty(LEV)
 tempp               = diff(LEV);
 precision           = floor(log10(abs(tempp)));
 precision(tempp==0) = 0;
 precision           = min(precision)-tol;
 % Rounds:
 LEV(1)       = floor(LEV(1)      *10^(-precision))*10^precision;
 LEV(2:end-1) = round(LEV(2:end-1)*10^(-precision))*10^precision;
 LEV(end)     = ceil(LEV(end)     *10^(-precision))*10^precision;
elseif ~isempty(Jlev)
 tempp               = diff(Jlev);
 if isempty(tempp)
  tempp              = Jlev;
 end
 precision           = floor(log10(abs(tempp)));
 precision(tempp==0) = 0;
 precision           = min(precision)-tol;
 % Rounds:
 if length(Jlev)==1
  Jlev          = round(Jlev*10^(-precision))*10^precision;
 else
  Jlev(1)       = floor(Jlev(1)      *10^(-precision))*10^precision;
  Jlev(2:end-1) = round(Jlev(2:end-1)*10^(-precision))*10^precision;
  Jlev(end)     = ceil(Jlev(end)     *10^(-precision))*10^precision;
 end
else
 tempp               = CLIM;
 precision           = floor(log10(abs(tempp)));
 precision(tempp==0) = 0;
 precision           = min(precision)-tol;
end
% Rounds:
CLIM(1) = floor(CLIM(1)*10^(-precision))*10^precision;
CLIM(2) =  ceil(CLIM(2)*10^(-precision))*10^precision;
warning(temp.state,'MATLAB:log:logOfZero')

% Completes levels when only join levels are specified:
if ~isempty(Jlev)
 cedge = CLIM;
 % First limit:
 if cedge(1)<=Jlev(1)
  if ~isempty(WID)
   cedge(1) = Jlev(1);
   if WID(1)~=0
    cedge(1) = cedge(1) - WID(1)*ceil((Jlev(1)-CLIM(1))/WID(1));
   end
  else
   % continue
  end
 else
  if (Ncmaps==2)
   cedge(1) = Jlev(1);
  else
   for k = 2:length(Jlev)
    if cedge(1)<=Jlev(k)
     cedge(1) = Jlev(k-1);
     break
    else
     Jlev(k-1) = Jlev(k);
    end
   end
  end
 end
 % Last limit:
 if cedge(2)>=Jlev(end)
  if ~isempty(WID)
   cedge(2) = Jlev(end);
   if WID(end)~=0
    cedge(2) = cedge(2) + WID(end)*ceil((CLIM(2)-Jlev(end))/WID(end));
   end
  else
   % continue
  end
 else
  if (Ncmaps==2)
   cedge(2) = Jlev(end);
  else
   for k = length(Jlev)-1:-1:1
    if cedge(2)>=Jlev(k)
     cedge(2) = Jlev(k+1);
     break
    else
     Jlev(k+1) = Jlev(k);
    end
   end
  end
 end
 % New Levels:
 LEV = [cedge(1); Jlev; cedge(2)];
 
end

% Gets colorband width and sets WID:
if ~isempty(Ncol)
 if isempty(WID)
  % Treats all colorbands with equal widths:
  Cwid = diff(CLIM)/sum(abs(Ncol));
  Cwid = round(Cwid*10^(-(precision-1)))*10^(precision-1);
  WID  = repmat(Cwid,Ncmaps,1);
  LEV  = [CLIM(1); CLIM(1)+cumsum(abs(Ncol))*Cwid];
 else
  % Treats WID as colorbands withs relations:
  WID   = WID/min(WID(WID~=0));
  Ncol2 = WID.*Ncol;
  Cwid  = diff(CLIM)/sum(abs(Ncol2));
  Cwid  = round(Cwid*10^(-(precision-1)))*10^(precision-1);
  WID   = WID*Cwid;
  LEV   = [CLIM(1); CLIM(1)+cumsum(abs(Ncol2))*Cwid];
 end
elseif ~isempty(WID)
 % Gets colorband width:
 Cwid  = WID(1)*10^(-precision);
 for k = 2:Ncmaps
  Cwid = gcd(Cwid,WID(k)*10^(-precision));
 end
 Cwid  = Cwid*10^precision;
else
 % Gets relation between colomaps width:
 if isempty(LEV)
  r    = ones(Ncmaps,1);
  d    = diff(CLIM);
 else
  r         = diff(LEV);
  temp      = warning('off','MATLAB:log:logOfZero');
  precision = floor(log10(abs(r))); % r = Str.XXX x 10^precision.
  precision(r==0) = 0; % precision=0 if Ncol=0.
  warning(temp.state,'MATLAB:log:logOfZero')
  precision = min(precision)-tol;
  r  = round(r*10^(-precision));
  rgcd  = r(1);
  for k = 2:Ncmaps
   rgcd = gcd(rgcd,r(k));
  end
  r = r/rgcd;
  d = (LEV(end)-LEV(1));
 end
 % Gets colorband width:
 r    = r*ceil(Tcol/sum(r));
 Cwid = d/sum(r);
 WID  = repmat(Cwid,Ncmaps,1);
end

% Sets LEV when empty:
if isempty(LEV)
 LEV = linspace(CLIM(1),CLIM(2),Ncmaps+1)';
end

% Gets number of colors for each colormap:
Ncol2 = round(diff(LEV)/Cwid);
if ~isempty(Ncol)
 % continue
else
 Ncol = round(diff(LEV)./WID);
 Ncol(~isfinite(Ncol)) = 0;
 if ~all(Ncol==round(Ncol))
  error('CVARGAS:cmjoin:incorrectWidColor',...
   'Colorband do not match each colormap width. Modify LEV or WID.')
 end
end

% Generates the colormaps:
CMAP  = zeros(sum(abs(Ncol2)),3);
xband = zeros(sum(abs(Ncol2))+1,1);
tempr = [];
for k = 1:Ncmaps
 if Ncol(k)
  r          = sum(abs(Ncol2(1:k-1)))+(1:abs(Ncol2(k)));
  if Ncol(k)~=Ncol2(k)
   CMAP(r,:) = cmapping(Ncol2(k),cmapping(Ncol(k),CMAPS{k}),'discrete');
  else
   CMAP(r,:) = cmapping(Ncol(k),CMAPS{k});
  end
  tempr      = linspace(LEV(k),LEV(k+1),abs(Ncol2(k))+1)';
  xband(r)   = tempr(1:end-1); 
 end
end
if ~isempty(tempr)
 xband(end) = tempr(end);
end

% Cuts edges:
ind = find((xband>=CLIM(1)) & (xband<=CLIM(2)));
if (ind(1)~=1) && ~(any(xband==CLIM(1)))
 ind = [ind(1)-1; ind];
end
if (ind(end)~=length(ind)) && ~(any(xband==CLIM(2)))
 ind = [ind; ind(end)+1];
end
CMAP  = CMAP(ind(1:end-1),:);
clim2 = xband(ind([1 end]));


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

if ~nargout
 colormap(AX{:},CMAP)
 caxis(AX{:},clim2(:)');
 clear CMAP
else
 if ~isempty(AX)
  colormap(AX{:},CMAP)
  caxis(AX{:},clim2(:)');
 end
 CLIM = clim2;
 WID  = diff(LEV)./max([Ncol ones(Ncmaps,1)],[],2);
end


% [EOF]   cmjoin.m
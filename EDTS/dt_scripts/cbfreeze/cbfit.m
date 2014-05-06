function CBH = cbfit(varargin)
%CBFIT   Draws a colorbar with specific color bands between its ticks.
% 
%   SYNTAX:
%           cbfit
%           cbfit(NBANDS)               % May be LBANDS instead of NBANDS
%           cbfit(NBANDS,CENTER)
%           cbfit(...,MODE)
%           cbfit(...,OPT)
%           cbfit(CBH,...)
%     CBH = cbfit(...);
%
%   INPUT:
%     NBANDS - Draws a colorbar with NBANDS bands colors between each tick
%      or      mark or a colorband between the specifies level bands
%     LBANDS   (LBANDS=NBANDS).
%              DEFAULT: 5     
%     CENTER - Center the colormap to this CENTER reference.
%              DEFAULT: [] (do not centers)
%     MODE   - Specifies the ticks mode (should be before AP,AV). One of:
%                'manual' - Forces color ticks on the new bands. 
%                'auto'   - Do not forces
%              DEFAULT: 'auto'
%     OPT    - Normal optional arguments of the COLORBAR function (should
%              be the last arguments).
%              DEFAULT: none.
%     CBH    - Uses this colorbar handle instead of current one.
%
%   OUTPUT (all optional):
%     CBH  - Returns the colorbar axes handle.
%
%   DESCRIPTION:
%     Draws a colorbar with specified number of color bands between its
%     ticks by modifying the current colormap and caxis.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * Sets the color limits, CAXIS, and color map, COLORMAP, before using
%       this function. Use them after this function to get the
%       modifications.
%
%   EXAMPLE:
%     figure, surf(peaks+2), colormap(jet(14)), colorbar
%      title('Normal colorbar.m')
%     figure, surf(peaks+2),                    cbfit(2,0)
%      title('Fitted 2 color bands and centered on zero')
%     figure, surf(peaks+2), caxis([0 10]),     cbfit(4,8)
%      title('Fitted 4 color bands and centered at 8')
%
%   SEE ALSO:
%     COLORBAR
%     and
%     CBFREEZE, CMFIT by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbfit.m
%   VERSION: 2.1 (Sep 30, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released as COLORBARFIT.M. (Mar 11, 2008)
%   1.1      Fixed bug when CAXIS is used before this function. (Jul 01,
%            2008)
%   1.2      Works properly when CAXIS is used before this function. Bug
%            fixed on subfunction and rewritten code. (Aug 21, 2008)
%   2.0      Rewritten code. Instead of the COLORBAND subfunction, now uses
%            the CMFIT function. Changed its name from COLORBARFIT to
%            CBFIT. (Jun 08, 2008)
%   2.1      Fixed bug and help with CBH input. (Sep 30, 2009)

%   DISCLAIMER:
%   cbfit.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2008,2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Sets defaults:
NBANDS = 5;
CENTER = [];
MODE   = 'auto';            
CBH    = [];
pax    = [];        % Peer axes

% Checks if first argument is a handle: Fixed bug Sep 2009
if (~isempty(varargin) && (length(varargin{1})==1) && ...
  ishandle(varargin{1})) && strcmp(get(varargin{1},'Type'),'axes')
 if strcmp(get(varargin{1},'Tag'),'Colorbar')
  CBH = varargin{1};
 else
  warning('CVARGAS:cbfit:incorrectHInput',...
   'Unrecognized first input handle.')
 end
 varargin(1) = [];
end
 
% Reads NBANDS and CENTER:
if ~isempty(varargin) && isnumeric(varargin{1})
 if ~isempty(varargin{1})
  NBANDS = varargin{1};
 end
 if (length(varargin)>1) && isnumeric(varargin{2})
  CENTER = varargin{2};
  varargin(2) = [];
 end
 varargin(1) = [];
end

% Reads MODE:
if (~isempty(varargin) && (rem(length(varargin),2)==1))
 if (~isempty(varargin{1}) && ischar(varargin{1}))
  switch lower(varargin{1})
   case 'auto'  , MODE = 'auto';
   case 'manual', MODE = 'manual';
   otherwise % 'off', 'hide' and 'delete'
    warning('CVARGAS:cbfit:incorrectStringInput',...
     'No pair string input must be one of ''auto'' or ''manual''.')
  end
 end
 varargin(1) = [];
end

% Reads peer axes:
for k = 1:2:length(varargin)
 if ~isempty(varargin{k})
  switch lower(varargin{k})
   case 'peer', pax = varargin{k+1}; break
  end
 end
end
if isempty(pax)
 pax = gca;
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Generates a preliminary colorbar:
if isempty(CBH)
 CBH = colorbar(varargin{:});
end

% Gets limits and orientation:
s     = 'Y';
ticks = get(CBH,[s 'Tick']);
if isempty(ticks)             
 s     = 'X';
 ticks = get(CBH,[s 'Tick']);
end
zlim = get(CBH,[s 'Lim']);

% Gets width and ref:
if ~isempty(NBANDS)

 NL = length(NBANDS);
 
 if (NL==1)
  
  % Force positive integers:
  NBANDS = round(abs(NBANDS));
 
  % Ignores ticks outside the limits:
  if zlim(1)>ticks(1)
   ticks(1) = [];
  end
  if zlim(2)<ticks(end)
   ticks(end) = [];
  end

  % Get the ticks step and colorband:
  tstep = ticks(2)-ticks(1);
  WIDTH = tstep/NBANDS;
  
  % Sets color limits
  if strcmp(get(pax,'CLimMode'),'auto')
   caxis(zlim);
  end
  
  % Forces old colorbar ticks: 
  set(CBH,[s 'Lim'],zlim,[s 'Tick'],ticks)
  
  % Levels:
  if strcmp(MODE,'manual')
   LBANDS = [fliplr(ticks(1)-WIDTH:-WIDTH:zlim(1)) ticks(1):WIDTH:zlim(2)];
  end
  
 else
  
  % Nonlinear colorbar:
  ticks = NBANDS;
  WIDTH = ticks;
  
  % Scales to CLIM:
  if strcmp(get(pax,'CLimMode'),'manual')
   ticks = ticks-ticks(1);
   ticks = ticks/ticks(end);
   ticks = ticks*diff(zlim) + zlim(1);
  end
  zlim = [ticks(1) ticks(end)];
  caxis(pax,zlim)
  CBIH = get(CBH,'Children');
  
  % Change ticks:
  set(CBIH,[s 'Data'],ticks)
  
  % Sets limits:
  set(CBH,[s 'Lim'],zlim)
  
  % Levels:
  if strcmp(MODE,'manual')
   LBANDS = NBANDS;
  end
  
 end
 
 % Get reference mark
 if ~isempty(CENTER)
  REF    = CENTER;
  CENTER = true;
 else
  REF    = ticks(1);
  CENTER = false;
 end
  
end

% Fits the colormap and limits:
cmfit(get(get(pax,'Parent'),'Colormap'),zlim,WIDTH,REF,CENTER)

% Sets ticks:
if strcmp(MODE,'manual')
 set(CBH,[s 'Tick'],LBANDS)
end

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

if ~nargout
 clear CBH
end


% [EOF]   cbfit.m
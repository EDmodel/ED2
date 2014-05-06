function CBLH = cblabel(varargin)
%CBLABEL   Adds a label to the colorbar.
%
%   SYNTAX:
%            cblabel(LABEL)
%            cblabel(LABEL,..,TP,TV);
%            cblabel(H,...)
%     CBLH = cblabel(...);
%
%   INPUT:
%     LABEL - String (or cell of strings) specifying the colorbar label.
%     TP,TV - Optional text property/property value arguments (in pairs).
%             DEFAULT:  (none)
%     H     - Color bar or peer axes (see COLORBAR) or figure handle(s) to
%             search for a single color bar handle.
%             DEFAULT: gca (current axes color bar)
%
%   OUTPUT (all optional):
%     CBLH  - Returns the colorbar label handle(s).
%           - Labels modified on the colorbar of the current figure or
%             the one(s) specified by CBH.
%
%   DESCRIPTION:
%     This function sets the label of the colorbar(s) in the current
%     figure.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%
%   EXAMPLE:
%     figure, colorbar, cblabel(['           T, °C'],'Rotation',0)
%     figure
%      subplot(211), h1 = colorbar;
%      subplot(212), h2 = colorbar('Location','south');
%      cblabel([h1 h2],{'$1-\alpha$','$\beta^3$'},'Interpreter','latex')   
%
%   SEE ALSO: 
%     COLORBAR
%     and 
%     CBUNITS, CBHANDLE, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cblabel.m
%   VERSION: 2.0 (Jun 08, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Aug 21, 2008)
%   2.0      Minor changes. Added CBHANDLE dependency. (Jun 08, 2009)

%   DISCLAIMER:
%   cblabel.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2008,2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
cbappname = 'Frozen';         % Colorbar application data with fields:
                              % 'Location' from colorbar
                              % 'Position' from peer axes befor colorbar
                              % 'pax'      handle from peer axes.

% Sets defaults:
H     = get(get(0,'CurrentFigure'),'CurrentAxes');
LABEL = '';
TOPT  = {};
CBLH  = [];

% Number of inputs:
if nargin<1
 error('CVARGAS:cblabel:incorrectNumberOfInputs',...
        'At least one input is required.')
end

% Looks for H:
if nargin && ~isempty(varargin{1}) && all(ishandle(varargin{1}))
 H = varargin{1};
 varargin(1) = [];
end

% Looks for CBH:
CBH = cbhandle(H);
if isempty(CBH), if ~nargout, clear CBLH, end, return, end

% Looks for LABEL:
if ~isempty(varargin) && (ischar(varargin{1}) || iscellstr(varargin{1}))  
 LABEL = varargin{1};
 varargin(1) = [];
end

% Forces cell of strings:
if ischar(LABEL)
 % Same label to all the color bars:
 LABEL = repmat({LABEL},length(CBH),1);
elseif iscellstr(LABEL) && (length(LABEL)==length(CBH))
  % Continue...
else
 error('CVARGAS:cblabel:incorrectInputLabel',...
        ['LABEL must be a string or cell of strings of equal size as ' ...
         'the color bar handles: ' int2str(length(CBH)) '.'])
end

% OPTIONAL arguments:
if ~isempty(varargin)
 TOPT = varargin;
end
if length(TOPT)==1
 TOPT = repmat({TOPT},size(CBH));
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
% NOTE: Only CBH, LABEL and TOPT are needed.

% Applies to each colorbar:
CBLH = repmat(NaN,size(CBH));
for icb = 1:length(CBH)
 
 % Searches for label location:
 try 
  % Normal colorbar:
  location = get(CBH(icb),'Location');
 catch
  % Frozen colorbar:
  location = getappdata(CBH(icb),cbappname);
  location = location.Location;
 end
 switch location(1)
  case 'E', as  = 'Y';
  case 'W', as  = 'Y';
  case 'N', as  = 'X';
  case 'S', as  = 'X';
 end
 % Gets label handle:
 CBLH(icb) = get(CBH(icb),[as 'Label']);
 % Updates label:
 set(CBLH(icb),'String',LABEL{icb},TOPT{:});
 
end

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Sets output:
if ~nargout
 clear CBLH
end


% [EOF]   cblabel.m
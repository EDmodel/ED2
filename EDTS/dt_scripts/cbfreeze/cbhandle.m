function CBH = cbhandle(varargin)
%CBHANDLE   Handle of current colorbar axes.
%
%   SYNTAX:
%     CBH = cbhandle;
%     CBH = cbhandle(H);
%
%   INPUT:
%     H - Handles axes, figures or uipanels to look for colorbars.
%         DEFAULT: gca (current axes)
%
%   OUTPUT:
%     CBH - Color bar handle(s).
%
%   DESCRIPTION:
%     By default, color bars are hidden objects. This function searches for
%     them by its 'axes' type and 'Colorbar' tag.
%    
%   SEE ALSO:
%     COLORBAR
%     and
%     CBUNITS, CBLABEL, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbhandle.m
%   VERSION: 1.1 (Aug 20, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.1      Fixed bug with colorbar handle input. (Aug 20, 2009)

%   DISCLAIMER:
%   cbhandle.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
axappname = 'FrozenColorbar'; % Peer axes application data with frozen
                              % colorbar handle.

% Sets default:
H = get(get(0,'CurrentFigure'),'CurrentAxes');

if nargin && ~isempty(varargin{1}) && all(ishandle(varargin{1}))
 H = varargin{1};
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Looks for CBH:
CBH = [];
% set(0,'ShowHiddenHandles','on')
for k = 1:length(H)
 switch get(H(k),'type')
  case {'figure','uipanel'}
   % Parents axes?:
   CBH = [CBH; ...
    findobj(H(k),'-depth',1,'Tag','Colorbar','-and','Type','axes')];
  case 'axes'
   % Peer axes?:
   hin  = double(getappdata(H(k),'LegendColorbarInnerList'));
   hout = double(getappdata(H(k),'LegendColorbarOuterList'));
   if     (~isempty(hin)  && ishandle(hin))
    CBH = [CBH; hin];
   elseif (~isempty(hout) && ishandle(hout))
    CBH = [CBH; hout];
   elseif isappdata(H(k),axappname)
    % Peer from frozen axes?:
    CBH = [CBH; double(getappdata(H(k),axappname))];
   elseif strcmp(get(H(k),'Tag'),'Colorbar') % Fixed BUG Aug 2009
    % Colorbar axes?
    CBH = [CBH; H(k)];
   end
  otherwise
   % continue
 end
end
% set(0,'ShowHiddenHandles','off')


% [EOF]   cbhandle.m
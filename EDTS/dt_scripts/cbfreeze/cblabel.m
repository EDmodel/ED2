function CBLH = cblabel(varargin)
%CBLABEL   Adds a label to the colorbar.
%
%   SYNTAX:
%            cblabel(LABEL)
%            cblabel(LABEL,..,TP,TV);
%            cblabel(H,...)
%     CBLH = cblabel(...);
%
%   INPUTS:
%     LABEL - String (or cell of strings) specifying the colorbar label.
%     TP,TV - Optional text property/property value arguments (in pairs).
%             DEFAULT:  (none)
%     H     - Color bar or peer axes (see COLORBAR) or figure handle(s) to
%             search for a single color bar handle.
%             DEFAULT: gca (current axes color bar)
%
%   OUTPUTS (all optional):
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
%   VERSION: 2.2 (Jul 03, 2014) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>)
%   MATLAB:  8.2.0.701 (R2013b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com
%   REVISIONS:
%   1.0      Released. (Aug 21, 2008)
%   2.0      Minor changes. Added CBHANDLE dependency. (Jun 08, 2009)
%   2.1      Changed application data to 'cbfreeze'. (Jun 05, 2014)
%   2.2      Fixed small bug with input reading. (Jul 03, 2014)
%   DISCLAIMER:
%   cblabel.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.
%   Copyright (c) 2008-2014 Carlos Adrian Vargas Aguilera
% INPUTS CHECK-IN
% -------------------------------------------------------------------------
% Parameters:
appName = 'cbfreeze';
% Sets defaults:
H     = get(get(0,'CurrentFigure'),'CurrentAxes');
LABEL = '';
TOPT  = {};
% Number of inputs:
assert(nargout<=1, 'CVARGAS:cblabel:tooManyOutputs',...
    'At most 1 output is allowed.')
% Looks for H: Version 2.2
if ~isempty(varargin) && ~isempty(varargin{1}) && ...
        all(reshape(ishandle(varargin{1}),[],1))
    H = varargin{1};
    varargin(1) = [];
end
% Looks for CBH:
CBH  = cbhandle(H,'force');
Ncbh = length(CBH);
% Looks for LABEL:
if ~isempty(varargin) && (ischar(varargin{1}) || iscellstr(varargin{1}))
    LABEL = varargin{1};
    varargin(1) = [];
end
% Forces cell of strings:
if ischar(LABEL)
    % Same label to all the colorbars:
    LABEL = repmat({LABEL},Ncbh,1);
elseif iscellstr(LABEL) && (length(LABEL)==Ncbh)
    % Continue...
else
    error('CVARGAS:cblabel:incorrectInputLabel',...
        ['LABEL must be a string or cell of strings of equal size as ' ...
        'the color bar handles: %d.'],Ncbh)
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
% Applies to each colorbar:
CBLH = NaN(1,Ncbh);
for icb = 1:Ncbh
    
    % Colorbar handle
    cbh = double(CBH(icb));
    
    % Searches for label location:
    try
        % Normal colorbar:
        location = get(cbh,'Location');
    catch
        % Frozen colorbar:
        location = getappdata(cbh,appName);
        location = location.cbLocation;
    end
    switch location(1)
        case 'E', XYstr  = 'Y';
        case 'W', XYstr  = 'Y';
        case 'N', XYstr  = 'X';
        case 'S', XYstr  = 'X';
    end
    
    % Gets label handle:
    CBLH(icb) = get(cbh,[XYstr 'Label']);
    
    % Updates label:
    set(CBLH(icb),'String',LABEL{icb},TOPT{:});
    
end
% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------
% Sets output:
if ~nargout
    clear CBLH
else
    CBLH(~ishandle(CBLH)) = [];
end
end
% [EOF] CBLABEL.M by Carlos A. Vargas A.

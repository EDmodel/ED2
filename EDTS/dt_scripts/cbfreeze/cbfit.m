function CBH = cbfit(varargin)
%CBFIT   Changes COLORMAP and CAXIS to fit between colorbar's ticks.
%
%   SYNTAX:
%           cbfit
%           cbfit(NBANDS)               % Or LBANDS
%           cbfit(NBANDS,CENTER)
%           cbfit(...,'force')
%           cbfit(H,...)
%     CBH = cbfit(...);
%
%   INPUTS:
%     NBANDS  - If scalar, fits NBANDS colors between each colorbar ticks.
%               If vector, fits a single color between specified LBANDS
%               ticks.
%               DEFAULT: 5
%     CENTER  - Centers the colormap to this reference.
%               DEFAULT: [] (do not centers)
%     'force' - By default, CBFIT changes the actual ticks. Use this
%               option to use the actual ones or the ones at LBANDS.
%     H       - Uses this colorbar handle (or peer axes, or figure)
%               instead of current one.
%
%   OUTPUTS (all optional):
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
%     Z = peaks; Z(Z<0) = Z(Z<0)*3;
%     figure(1)
%       surf(Z)
%       colorbar
%       title({'Normal "ugly" COLORBAR'; ...
%              'with colorbar ticks anywhere'; ...
%              'and "green" centered at -5!'})
%     figure(2)
%       surf(Z)
%       cbfit(4,0)
%       title({'CBFIT(4,0)'; ...
%              'Fitted 4 color bands between ticks'; ...
%              'and "green" centered at zero'})
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
%   VERSION: 3.0 (Jun 05, 2014) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>)
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
%   3.0      Rewritten code. Changed optional input 'manual' to 'force'.
%            Now accepts axes and figure handles as inputs. (Jun 05, 2014)
%   DISCLAIMER:
%   cbfit.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.
%   Copyright (c) 2008-2014 Carlos Adrian Vargas Aguilera
% INPUTS CHECK-IN
% -------------------------------------------------------------------------
% Sets defaults:
NBANDS = 5;
LBANDS = [];
CENTER = [];
FORCE  = false;
H      = get(get(0,'CurrentFigure'),'CurrentAxes');
% Checks inputs/outputs number:
assert(nargin<=4, 'CVARGAS:cbfit:tooManyInputs',...
    'At most 4 inputs are allowed.')
assert(nargout<=1, 'CVARGAS:cbfit:tooManyOutputs',...
    'At most 1 output is allowed.')
% Reads FORCE:
if ~isempty(varargin) && ~isempty(varargin{end}) && ...
        ischar(varargin{end})
    if numel(varargin{end})==size(varargin{end},2)
        switch lower(varargin{end})
            case {'force','forc','for','fo','f'}
                FORCE = true;
        end
    end
    varargin(end) = [];
end
% Reads NBANDS and CENTER:
if ~isempty(varargin) && isnumeric(varargin{end})
    if (length(varargin)>1) && isnumeric(varargin{end-1})
        CENTER = varargin{end};
        varargin(end) = [];
    end
    if numel(varargin{end})==1
        NBANDS = varargin{end};
        LBANDS = [];
    else
        LBANDS = reshape(varargin{end},[],1);
        NBANDS = [];
    end
    varargin(end) = [];
end
% Reads H:
if (length(varargin)==1) && ~isempty(varargin{1}) && ...
        all(reshape(ishandle(varargin{1}),[],1))
    H = varargin{1};
end
% Gets colorbar handles or creates them:
CBH = cbhandle(H,'force');
Ncbh = length(CBH);
% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
for icb = 1:Ncbh
    % Generates a preliminary colorbar:
    
    % Colorbar handle:
    cbh = double(CBH(icb));
    
    % Gets limits and orientation:
    XYstr = 'Y';
    ticks = get(cbh,[XYstr 'Tick']);
    if isempty(ticks)
        XYstr = 'X';
        ticks = get(cbh,[XYstr 'Tick']);
    end
    XYLim = get(cbh,[XYstr 'Lim']);
    tempp = ticks;
    % Gets peer axes
    peer = cbhandle(cbh,'peer');
    
    % Gets width and ref:
    if ~isempty(NBANDS)
        
        % Force positive integers:
        NBANDS = round(abs(NBANDS));
        
        % Ignores ticks outside the limits:
        if XYLim(1)>ticks(1)
            ticks(1) = [];
        end
        if XYLim(2)<ticks(end)
            ticks(end) = [];
        end
        
        % Get the ticks step and colorband:
        tstep  = ticks(2)-ticks(1);
        WIDTH  = tstep/NBANDS;
        LBANDS = ticks;
        % LBANDS = [fliplr(ticks(1)-WIDTH:-WIDTH:XYLim(1)) ...
        %     ticks(1):WIDTH:XYLim(2)];
        
        % Sets color limits
        if strcmp(get(peer,'CLimMode'),'auto')
            caxis(XYLim);
        end
        
        % Forces old colorbar limits:
        set(cbh,[XYstr 'Lim'],XYLim) % ,[XYstr 'Tick'],ticks)
        
    else % ~isempty(LBANDS)
            
        % Nonlinear colorbar:
        ticks = LBANDS;
        WIDTH = ticks;
        
        % Scales to CLIM:
        if strcmp(get(peer,'CLimMode'),'manual')
            ticks = ticks-ticks(1);
            ticks = ticks/ticks(end);
            ticks = ticks*diff(XYLim) + XYLim(1);
        end
        XYLim = [ticks(1) ticks(end)];
        caxis(peer,XYLim)
        CBIH = get(cbh,'Children');
        
        % Change ticks:
        set(CBIH,[XYstr 'Data'],ticks)
        
        % Sets limits:
        set(cbh,[XYstr 'Lim'],XYLim)
        
    end
        
    % Get reference mark
    if ~isempty(CENTER)
        REF    = CENTER;
        CENTER = true;
    else
        REF    = ticks(1);
        CENTER = false;
    end
        
    % Fits the colormap and limits:
    fig  = get(peer,'Parent');
    CMAP = colormap(fig);
    cmfit(CMAP,XYLim,WIDTH,REF,CENTER);
            
    % Sets ticks:
    if FORCE
        set(cbh,[XYstr 'Tick'],LBANDS)
    end
end
% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------
if ~nargout
    clear CBH
end
end
% [EOF] CBFIT.M by Carlos A. Vargas A.

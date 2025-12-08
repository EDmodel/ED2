function CBH = cbunits(varargin)
%CBUNITS   Adds units (and ISU prefixes) to the colorbar ticklabels.
%
%   SYNTAX:
%           cbunits ºC             % Just an example for celsius degrees.
%           cbunits ºC prefix      % Uses ISU's prefix instead of 'x 10^y'.
%           cbunits ºC off         % Resets the labels.
%           cbunits(H,'°C',...)    % Uses H axes instead of current ones.
%     CBH = cbunits(...);
%
%   INPUTS:
%     UNITS - Units to be added in the current COLORBAR tick labels. May be
%             a string or a cell of strings if several colorbars are used.
%     H     - Use this colorbar(s) handle(s) (or peer axes (see COLORBAR)
%             or figure handle(s)) instead of the current one.
%             DEFAULT: gca (uses the current one or creates one)
%     '...' - String optional inputs as follows:
%              ------------------------------------------------------------
%               OPTION       DESCRIPTION
%              ------------------------------------------------------------
%               'prefix'     Uses International System Units' prefixes.
%                            Like 'M' (mega), 'd' (deci), 'p' (pico),
%                            instead of the small scientific notation text
%                            'x10^6', 'x10^-2', 'x10^-12', respectively.
% 
%               'off'        Undoes previous use of CBUNITS.
%              ------------------------------------------------------------
%              Other char options are listed on the NOTEs.
%
%   OUTPUTS (all optional):
%     CBH   - Returns the colorbar handle(s).
%             DEFAULT: Not returned if not required.
%           - Ticklabels modified on the colorbar of the current axes or
%             the one(s) specified by CBH.
%
%   DESCRIPTION:
%     This function adds units to the current colorbar, by writting them
%     on the highest ticklabel. It works also with logarithmic scales.
%
%     When a scientific notation is required, MATLAB normally put it in a
%     small textbox like 'x10^6'; this function will put it within the unit
%     as '°F/1e6' (dividing the unit), for example; or as 'M°F' if 'prefix'
%     was used. If no unit was provided ('') then uses '1e6' (multiplying
%     the label).
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * When more than one colorbar handle is given or founded and a single
%       UNITS string is given, it is applied to all of them.
%     * If not COLORBAR is found on the H handle(s), it creates one on each
%       of the axes.
%     * Use a cell of strings for UNITS when more than one colorbar handles
%       are given in order to give to each one their proper units. This
%       also works when the handles are founded but the units order is
%       confusing and not recommended.
%     * Once applied, CAXIS shouldn't be used.
%     * Instead of 'prefix' you may use a specific prefix like 'mili',
%       'mega', 'kilo', 'deci', 'nano', etc. But this functionality is not 
%       recommended because you can get huge numbers when used incorrectly.
%     * Always an extra space is put between the label and the unit, use
%       'nospace' extra input to avoid it.
%     * Use 'upE' extra input to use uppercase 'E' instead of 'e' for the
%       scientific notation. 
%     * Use 'normal' extra input to avoid forcing all labels with the same
%       decimals, that is, to use normal MATLAB behaviour.
%     * Give extra inputs as a list of chars, for example:
%         >> cbunits °C prefix nospace
%     * Use '1' for adimensional values, for example:
%         >> cbunits 1 prefix
%
%   EXAMPLE: (COLORMAP with MEGAWATTS units)
%       contour(peaks(30)*1e7)
%       cbunits W prefix
%
%   SEE ALSO:
%     COLORBAR
%     and
%     CBLABEL, CBHANDLE, CBFREEZE by Carlos Vargas
%     ticks http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbunits.m
%   VERSION: 4.4 (Jul 03, 2014) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>)
%   MATLAB:  8.2.0.701 (R2013b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com
%   REVISIONS:
%   1.0      Released. (Aug 21, 2008)
%   2.0      Minor changes. Added 'clear' option and CBHANDLE dependency.
%            (Jun 08, 2009)
%   3.0      Fixed bug when inserting units on lower tick and ticklabel
%            justification. Added SPACE option. (Sep 30, 2009)
%   4.0      Eliminated the 'clear' option. Changed SPACE logical option to
%            string 'nospace'. Included several options for scientiic
%            notation, like 'prefix' (for auto notation) or 'kilo', 'mega',
%            etc. (Jun 05, 2014)
%   4.1      If no unit is given, sets scientific notation like 'eN' 
%            instead of MATLAB's 'x10N'. Forces all labels to have the same
%            number of decimals. (Jun 14, 2014)  
%   4.2      Fixed several bugs with 'log' scale. (Jun 15, 2014)
%   4.3      Fixed small bug with '1' unit. (Jun 17, 2014)
%   4.4      Fixed small bug with 'e' char. Added extra options 'upE' and
%            'normal'. (Jul 03, 2014) 
%   DISCLAIMER:
%   cbunits.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.
%   Copyright (c) 2008-2014 Carlos Adrian Vargas Aguilera
% INPUTS CHECK-IN
% -------------------------------------------------------------------------
% Sets defaults:
UNITS   = '';
SPACE   = true;
PREFIX  = {'',[]};
OFF     = false;
EXP     = 'e';    % Letter for exponential notation: 'e' or 'E'.
FORCE   = 1==1;   % All labels with the same decimals (Version 4.1)
H       = get(get(0,'CurrentFigure'),'CurrentAxes');
% Constants:
mu      = char(181); % For helvetica FONT!
prepow  = [-24:3:-3 -2 -1 1 2 3:3:24];
prestr  = 'yzafpnmmcddhkMGTPEZY';
% Anonymous definitions:
numsubs = @(x,ind,val) ... % Numerical substitution, x(ind) = val
    subsasgn( x, struct( 'type', '()', 'subs', {{ind}}), val) + 0;
oom     = @(x) ...         % Order of magnitud of x
    numsubs( floor(log10( abs(x))),         x==0,   0) + ...
    numsubs(       zeros(size(x)) , ~isfinite(x), NaN);
sig     = @(x) ...         % Significant digits on x (Version 4.1)
    bsxfun( @(y,sigy) min( numsubs( sigy, sigy==0, Inf), ...
    cell2mat( regexp( arrayfun( @(sigz,z) sprintf( '%1.*e', ...
    min(sigz,15)-1, z), numsubs( sigy, isnan(sigy), 1), y, 'Uni', ...
    false), '((\d|\.)0*e[\+-][\d]{1,3})|(NaN|[-]?Inf)$')) ) - 1, ...
    abs(x), oom(x)-oom(eps(x)));
pre     = @(x) ...         % Precision of x
    oom(x) - sig(x) + 1;
sliceit = @(x,n) arrayfun(@(sind) x(1:sind),...
    regexp(x,['(?<=' x(1:n-1) '.*).']),'Uni',false);
% Checks inputs/outputs number:
assert(nargout<=1,'CVARGAS:cbunits:tooManyOutputs',...
    'At most 1 output is allowed.')
% Reads H:
if ~isempty(varargin) && ~isempty(varargin{1}) && ...
        all(reshape(ishandle(varargin{1}),[],1))
    H = varargin{1};
    varargin(1) = [];
end
% Reads UNITS:
if ~isempty(varargin) && (ischar(varargin{1}) || iscellstr(varargin{1}))
    UNITS = varargin{1};
    varargin(1) = [];
end
% Reads optional inputs:
while ~isempty(varargin)
    if ~ischar(varargin{1}) || (numel(varargin{1})~=size(varargin{1},2))
        varargin(1) = [];
        continue
    end
    switch lower(varargin{1})
        case sliceit('off',1)
            OFF    = true;
        case sliceit('nospace',3)
            SPACE  = false;
        case sliceit('normal',3)
            FORCE  = false;
        case sliceit('upe',1)
            EXP    = 'E';
        case sliceit('prefix',2)
            PREFIX = {'prefix',[]};
        case sliceit('yotta',3)
            PREFIX = {'Y',1e24};
        case sliceit('zetta',3)
            PREFIX = {'Z',1e21};
        case sliceit('exa',1)
            PREFIX = {'E',1e18};
        case sliceit('peta',2)
            PREFIX = {'P',1e15};
        case sliceit('tera',1)
            PREFIX = {'T',1e12};
        case sliceit('giga',1)
            PREFIX = {'G',1e9};
        case sliceit('mega',2)
            PREFIX = {'M',1e6};
        case sliceit('kilo',1)
            PREFIX = {'k',1e3};
        case sliceit('hecto',1)
            PREFIX = {'h',1e2};
        case sliceit('deca',4)
            PREFIX = {'da',1e1};
        case sliceit('deci',4)
            PREFIX = {'d',1e-1};
        case sliceit('centi',1)
            PREFIX = {'c',1e-2};
        case sliceit('mili',3)
            PREFIX = {'m',1e-3};
        case sliceit('micro',3)
            PREFIX = {mu,1e-6};
        case sliceit('nano',2)
            PREFIX = {'n',1e-9};
        case sliceit('pico',2)
            PREFIX = {'p',1e-12};
        case sliceit('femto',1)
            PREFIX = {'f',1e-15};
        case sliceit('atto',1)
            PREFIX = {'a',1e-18};
        case sliceit('zepto',3)
            PREFIX = {'z',1e-21};
        case sliceit('yocto',3)
            PREFIX = {'y',1e-24};
        otherwise
            warning('CAVARGAS:cbunits:UnrecognizedOptionalInput', ...
                'Optional input ''%s'' wasn''t recognized',varargin{1})
    end
    varargin(1) = [];
end
% Gets colorbar handles or creates them:
CBH  = cbhandle(H,'force');
Ncbh = length(CBH);
% Forces UNITS as a cell of strings:
if ischar(UNITS)
    if numel(UNITS)~=size(UNITS,2)
        error('CVARGAS:cbunits:IncorrectUnitsString',...
            'UNITS string must be a row vector.')
    end
    % Same units for all the colorbars:
    UNITS = repmat({UNITS},Ncbh,1);
elseif iscellstr(UNITS) && (length(UNITS)==Ncbh)
    % Continue...
else
    error('CVARGAS:cbunits:IncorrectInputUnits',...
        ['UNITS must be a string or cell of strings of equal size as ' ...
        'the colorbar handles: %d.'],Ncbh)
end
% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
% Applies to each colorbar:
for icb = 1:Ncbh
    
    cbh    = double(CBH(icb));
    units  = UNITS{icb};
    prefix = '';
    sufix  = '';
    space  = repmat(' ',1,SPACE);
    
    % Gets colorbar position and ticks:
    XYstr  = 'Y';
    ticks  = get(cbh,[XYstr 'Tick']);
    if isempty(ticks)
        XYstr = 'X';
        ticks = get(cbh,[XYstr 'Tick']);
    end
    
    % Deletes current labels:
    set(cbh,[XYstr 'TickLabelMode'],'auto')
    if OFF, continue, end
    
    % Gets scale:
    islog = strcmp(get(cbh,[XYstr 'Scale']),'log');
    
    % Gets current labels:
    labels = get(cbh,[XYstr 'TickLabel']);
    
    % Checks for scientific notation:
    switch PREFIX{1,1}
        
        case {'','prefix'}
            % Gets the order of magnitud:
            if ~islog
                ind = find(ticks~=0,1);
                ord = oom(ticks(ind)/str2double(labels(ind,:)));
                if ord~=0 % Fixed bug 4.0 version
                    if isempty(PREFIX{1,1})
                        
                        if isempty(units)
                            % Adds the scientific notation to the unit like
                            % 'e...': (Version 4.1)
                            if strcmp(units,'1') % (Version 4.3)
                                units = '';
                            end
                            prefix = sprintf('%s%d',EXP,ord); 
                        else
                            % Adds the scientific notation to the unit like
                            % 'UNIT/1e...': 
                            sufix = sprintf('/1%s%d',EXP,ord);
                        end
                        labels = ticks(:)/10^ord;
                    else
                        
                        % Adds the scientific notation as IS prefix like
                        % 'M...':
                        [~,ind] = min(abs(prepow-ord));
                        prenum = prepow(ind);
                        switch prenum
                            case 1
                                prefix = 'da';
                            case -6
                                prefix = mu;
                            otherwise
                                prefix = prestr(ind);
                        end
                        labels = ticks(:)/10^prenum;
                        
                    end
                end
            elseif strcmp(units,'1') % (Version 4.2)
                units = ''; 
            end
            
        otherwise
            % Uses specific IS prefix:
            if strcmp(units,'1') % (Version 4.1)
                units = ''; 
            end
            prefix = PREFIX{1,1};
            labels = ticks(:)/PREFIX{1,2};   
            
    end
    
    % Changes labels from numeric to char: (Version 4.1)
    if FORCE
        % All labels with the same decimals:
        if ischar(labels)
            labels = ticks(:);
        end
        fmt = '';
        if islog % (Version 4.2)
            labels = log10(labels);
            fmt = ['1' EXP]; % (Version 4.4)
        end
        ndec   = abs(min(min(pre(labels(labels~=0))),0));
        fmt    = sprintf('%s%%.%df',fmt,ndec);
        labels = num2str(labels,fmt);
    elseif ~ischar(labels) || islog
        fmt = '';
        if islog
            labels = log10(ticks);
            fmt = ['1' EXP]; % (Version 4.4)
        end
        % Integers without decimals:
        ndec   = abs(min(pre(labels),0.*labels));
        fmt    = [sprintf('%s%%.',fmt), ...
            strjoin( arrayfun( @(x) sprintf( '%d', x), ...
            ndec, 'Uni', false)', ['f\\n' fmt '%.']), 'f'];
        labels = num2str(labels,fmt);
    end
    
    % Looks where to put the units, at the first label or the last one:
    Nt  = size(labels,1);
    loc = Nt; % Fixed bug 3.0 version
    if strcmp(get(cbh,[XYstr 'Dir']),'reverse') % (Version 4.1)
        loc = 1;
    end
    
    % Adds prefix units and sufix with the units:
    Nu = length(units);
    Ne = length(space);
    Ns = length(sufix);
    Np = length(prefix);
    labels = [strjust(labels,'right') repmat(' ',Nt,Ne+Np+Nu+Ns)];
    labels(loc,end-Ne-Np-Nu-Ns+1:end) = [space prefix units sufix];
    
    % Gets labels' justification:
    JUST = 'center';
    if strcmp(XYstr,'Y') % Fixed bug 3.0 version
        if strcmp(get(cbh,[XYstr 'AxisLocation']),'right')
            JUST = 'left';
        else
            JUST = 'right';
        end
    end
    
    % Sets labels on the colormap:
    set(cbh,[XYstr 'TickLabel'],strjust(labels,JUST))
    
end % MAIN LOOP
% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------
% Sets output:
if ~nargout
    clear CBH
end
end
% [EOF] CBUNITS.M by Carlos A. Vargas A.

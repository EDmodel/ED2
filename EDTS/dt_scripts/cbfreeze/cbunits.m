function CBH = cbunits(varargin)
%CBUNITS   Adds units to the colorbar ticklabels.
%
%   SYNTAX:
%           cbunits(UNITS)
%           cbunits(UNITS,SPACE)
%           cbunits -clear
%           cbunits(H,...)
%     CBH = cbunits(...);
%
%   INPUT:
%     UNITS - String (or cell of strings) with the colorbar(s) units or
%             '-clear' to eliminate any unit. 
%     SPACE - Logical indicating whether an space should be put between
%             quantity and units. Useful when using '3°C', for example.
%             DEFAULT: true (use an space)
%     H     - Colorbar, or peer axes (see COLORBAR) or figure handle(s) to
%             search for color bars. 
%             DEFAULT: gca (current axes color bar)
%
%   OUTPUT (all optional):
%     CBH   - Returns the colorbar handle(s).
%             DEFAULT: Not returned if not required.
%           - Ticklabels modified on the colorbar of the current axes or
%             the one(s) specified by CBH.
%
%   DESCRIPTION:
%     This function adds units to the current colorbar, by writting them
%     after the biggest ticklabel.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * Scientific notation is included in the units (if any).
%     * When more than one colorbar handle is given or founded and a single
%       UNITS string is given, it is applied to all of them.
%     * Use a cell of strings for UNITS when more than one colorbar handles
%       are given in order to give to each one their proper units. This
%       also works when the handlesare founded but the units order is
%       confusing and not recommended.
%     * Once applied, CAXIS shouldn't be used.
%     * To undo sets the ticklabelmode to 'auto'.
%
%   EXAMPLE:
%     % Easy to use:
%       figure, caxis([1e2 1e8]), colorbar, cbunits('°F',false)
%     % Vectorized:
%       figure
%       subplot(211), h1 = colorbar;
%       subplot(212), h2 = colorbar;
%       cbunits([h1;h2],{'°C','dollars'},[false true])
%     % Handle input:
%       figure
%       subplot(211), colorbar;
%       subplot(212), colorbar('Location','North');
%       caxis([1e2 1e8])
%       cbunits(gcf,'m/s')
%
%   SEE ALSO: 
%     COLORBAR
%     and
%     CBLABEL, CBHANDLE, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbunits.m
%   VERSION: 3.0 (Sep 30, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Aug 21, 2008)
%   2.0      Minor changes. Added 'clear' option and CBHANDLE dependency.
%            (Jun 08, 2009)
%   3.0      Fixed bug when inserting units on lower tick and ticklabel
%            justification. Added SPACE option. (Sep 30, 2009)

%   DISCLAIMER:
%   cbunits.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2008,2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Sets defaults:
H     = get(get(0,'CurrentFigure'),'CurrentAxes');
UNITS = '';
SPACE = true;

% Checks inputs/outputs number:
if     nargin<1
 error('CVARGAS:cbunits:notEnoughInputs',...
  'At least 1 input is required.')
elseif nargin>3
 error('CVARGAS:cbunits:tooManyInputs',...
  'At most 3 inputs are allowed.')
elseif nargout>1
 error('CVARGAS:cbunits:tooManyOutputs',...
  'At most 1 output is allowed.')
end

% Looks for H:
if nargin && ~isempty(varargin{1}) && all(ishandle(varargin{1}))
 H = varargin{1};
 varargin(1) = [];
end

% Looks for CBH:
CBH = cbhandle(H);
if isempty(CBH), if ~nargout, clear CBH, end, return, end

% Looks for UNITS:
if ~isempty(varargin) && ~isempty(varargin{1}) && ...
  (ischar(varargin{1}) || iscellstr(varargin{1}))  
 UNITS = varargin{1};
 varargin(1) = [];
end
if isempty(UNITS), if ~nargout, clear CBH, end, return, end

% Forces cell of strings:
if ischar(UNITS)
 if numel(UNITS)~=size(UNITS,2)
  error('CVARGAS:cbunits:IncorrectUnitsString',...
        'UNITS string must be a row vector.')
 end
 % Same units to all the color bars:
 UNITS = repmat({UNITS},length(CBH),1);
elseif iscellstr(UNITS) && (length(UNITS)==length(CBH))
  % Continue...
else
 error('CVARGAS:cbunits:IncorrectInputUnits',...
        ['UNITS must be a string or cell of strings of equal size as ' ...
         'the color bar handles: ' int2str(length(CBH)) '.'])
end

% Looks for SPACE:
Nunits = length(UNITS);
if ~isempty(varargin) && ~isempty(varargin{1}) && ...
  ((length(varargin{1})==1) || (length(varargin{1})==Nunits))  
 SPACE = varargin{1};
end
SPACE = logical(SPACE);

% Forces equal size of SPACE and UNITS.
if (length(SPACE)==1) && (Nunits~=1)
 SPACE = repmat(SPACE,Nunits,1);
end


% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
% Note: Only CBH and UNITS are required.

% Applies to each colorbar:
for icb = 1:length(CBH)
 
 units  = UNITS{icb};
 space  = SPACE(icb);
 cbh    = CBH(icb);
 append = [];
 
 % Gets tick labels:
 as  = 'Y';
 at  = get(cbh,[as 'Tick']);
 if isempty(at)
  as = 'X';
  at = get(cbh,[as 'Tick']);
 end
 
 % Checks for elimitation:
 if strcmpi(units,'-clear')
  set(cbh,[as 'TickLabelMode'],'auto')
  continue
 end

             set(cbh,[as 'TickLabelMode'],'manual');
 old_ticks = get(cbh,[as 'TickLabel']);

 % Adds scientific notation:
 if strcmp(get(cbh,[as 'Scale']),'linear')
  ind = 1;
  if at(ind)==0
   ind = 2;
  end
  o  = log10(abs(at(ind)/str2double(old_ticks(ind,:))));
  sg = '';
  if at(ind)<0, sg = '-'; end
  if o>0
   append = [' e' sg int2str(o) ''];
  end
 end
 
 % Updates ticklabels:
 Nu = length(units);
 Na = length(append);
 Nt = size(old_ticks,1);
 loc = Nt; % Fixed bug, Sep 2009
 if (strcmp(as,'Y') && ((abs(at(1))>abs(at(Nt))) && ...
    (length(fliplr(deblank(fliplr(old_ticks( 1,:))))) > ...
     length(fliplr(deblank(fliplr(old_ticks(Nt,:)))))))) || ...
     (strcmp(as,'X') && strcmp(get(cbh,[as 'Dir']),'reverse'))
  loc = 1; 
 end
 new_ticks  = [old_ticks repmat(' ',Nt,Nu+(Na-(Na>0))+space)];
 new_ticks(loc,end-Nu-Na-space+1:end) = [append repmat(' ',1,space) units];
 if strcmp(as,'Y') % Fixed bug, Sep 2009
  if strcmp(get(cbh,[as 'AxisLocation']),'right')
   new_ticks = strjust(new_ticks,'left');
  else
   new_ticks = strjust(new_ticks,'right');
  end
 else
  new_ticks = strjust(new_ticks,'center');
 end
 set(cbh,[as 'TickLabel'],new_ticks)
 
end % MAIN LOOP


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Sets output:
if ~nargout
 clear CBH
end


% [EOF]   cbunits.m
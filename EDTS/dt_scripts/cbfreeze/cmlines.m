function [HL,CLIN] = cmlines(varargin)
% CMLINES   Change the color of plotted lines using the colormap.
%
%   SYNTAX:
%                 cmlines
%                 cmlines(CMAP)
%                 cmlines(H,...)
%     [HL,CLIN] = cmlines(...);
%   
%   INPUT:
%     CMAP - Color map name or handle to be used, or a Nx3 matrix of colors
%            to be used for each of the N lines or color char specifiers.
%            DEFAULT: jet.
%     H    - Handles of lines or from a axes to search for lines or from
%            figures to search for exes. If used, must be the first input.
%            DEFAULT: gca (sets colors for lines in current axes)
%
%   OUTPUT (all optional):
%     HL   - Returns the handles of lines. Is a cell array if several axes
%            handle were used as input.
%     CLIN - Returns the RGB colors of the lines. Is a cell array if
%            several axes handle were used as input.
%
%   DESCRIPTION:
%     Ths function colored the specified lines with the spectrum of the
%     given colormap. Ideal for lines on the same axes which means increase
%     (or decrease) monotonically.
%
%   EXAMPLE:
%     plot(reshape((1:10).^2,2,5))
%     cmlines
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%    
%   SEE ALSO:
%     PLOT and COLORMAP.
%     and
%     CMAPPING
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cmlines.m
%   VERSION: 1.0 (Jun 08, 2009) (<a href="matlab:web(['www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do',char(63),'objectType',char(61),'author',char(38),'objectId=1093874'])">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)

%   DISCLAIMER:
%   cmlines.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera

% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Set defaults:
HL   = {};
Ha   = gca;
CMAP = colormap;

% Checks number of inputs:
if nargin>2
 error('CVARGAS:cmlines:tooManyInputs', ...
  'At most 2 inputs are allowed.')
end
if nargout>2
 error('CVARGAS:cmlines:tooManyOutputs', ...
  'At most 2 outputs are allowed.')
end

% Checks handles of lines, axes or figure inputs:
Hl = [];
if (nargin~=0) && ~isempty(varargin{1}) && all(ishandle(varargin{1}(:))) ...
 && ((length(varargin{1})>1) || ~isa(varargin{1},'function_handle'))
 Ha = [];
 for k = 1:length(varargin{1})
  switch get(varargin{1}(k),'Type')
   case 'line'
    Hl = [Hl varargin{1}(k)];
   case 'axes'
    Ha = [Ha varargin{1}(k)];
   case {'figure','uipanel'}
    Ha = [Ha findobj(varargin{1}(k),'-depth',1,'Type','axes',...
                      '-not',{'Tag','Colorbar','-or','Tag','legend'})];
   otherwise
     warning('CVARGAS:cmlines:unrecognizedHandleInput',...
      'Ignored handle input.')
  end
 end
 varargin(1) = [];
end

% Looks for CMAP input:
if nargin && ~isempty(varargin) && ~isempty(varargin{1})
 CMAP = varargin{1};
end

% Gets line handles:
if ~isempty(Hl)
 HL{1} = Hl;
end
if ~isempty(Ha)
 for k = 1:length(Ha)
  Hl = findobj(Ha(k),'Type','line');
  if ~isempty(Hl)
   HL{end+1} = Hl;
  end
 end
end
if isempty(HL)
 if ~nargout
  clear HL
 end
 return
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Sets color lines for each set of lines:
Nlines = length(HL);
CLIN   = cell(1,Nlines);
for k  = 1:length(HL)
 
 % Interpolates the color map:
 CLIN{k} = cmapping(length(HL{k}),CMAP);

 % Changes lines colors:
 set(HL{k},{'Color'},mat2cell(CLIN{k},ones(1,size(CLIN{k},1)),3))
 
end

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

if ~nargout
 clear HL
elseif Nlines==1
 HL   = HL{1};
 CLIN = CLIN{1};
end


% [EOF]   cmlines.m
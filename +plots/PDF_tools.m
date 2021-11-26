

function PDF_tools(mode,varargin)
% function plots.PDF_tools(mode, [do_PDF], ...)
% This is a combination of some PDF handling utilities, simplifying the
%   process of generating multi-page PDF output for data analyses. 
% 
% INPUTS 
% mode = 'setup', 'page', 'compile', or <array of figure handles>
% 
% Example usage :  
%  do_PDF = any(named('-pdf'));
%  plots.PDF_tools('setup',do_PDF);
% 
%  for ii = 1:num_pages
%     make_figure_page( ... )
%     plots.PDF_tools(gcf, 'page-%04d.ps', ii)
%  end
% 
%  plots.PDF_tools('compile',do_PDF,'my-analysed-data (%d).pdf')
%
% mode = 'page' and mode = <handle> do the same thing with the exception
%   that mode = <handle> has the side effect of calling pause(0.05) every
%   time plots.PDF_tools is called (regardless of whether the optional
%   do_PDF is true or false), and then closing <handle> after generating
%   the page (only if do_PDF was true). This was a very common code pattern
%   for PDF generation.
% 
% Calvin Eiber, Sander Pietersen, Paul Martin
% Last updated 2021 Nov 21

named = @(v) strncmpi(v,varargin,length(v));
get_  = @(v) varargin{find(named(v))+1};

% is do_PDF or the local equivalent (do_pdf, do.pdf, etc) false ?
if nargin > 1 && islogical(varargin{1}) 
  if any(ishandle(mode)), pause(0.05), end
  if ~varargin{1}, return, end
  varargin(1) = [];
end
 

persistent ps_folder
if isempty(ps_folder), ps_folder = 'make-pdf\'; end


switch lower(mode)
  case 'setup'
    
    if exist([tempdir ps_folder],'dir'), 
       rmdir([tempdir ps_folder],'s'); 
    end
    mkdir([tempdir ps_folder]);
    
  case 'page'
    
    if any(ishandle(varargin{1})), 
         fig_id = varargin{1}; varargin(1) = []; 
    else fig_id = gcf; 
    end
    
    nA = sum(varargin{1} == '%')+1; 
    if any(named('-n')), nA = get_('-n'); end
    if nA == numel(varargin), varargin{end+1} = '-move'; end
    
    page = [ps_folder sprintf(varargin{1:nA})]; 
    figs_to_ps(fig_id, page, varargin{nA+1:end}); 
    
  case {'compile','combine'}, 
    
    if nargin == 1,
         filename = tools.file('get','Results PDF (%d).pdf','next');
    elseif any(varargin{1} == '%')
         filename = tools.file('get',varargin{1},'next');
    else filename = tools.file(varargin{:});
    end
    
    combine_PStoPDF([tempdir ps_folder], filename), 
    if ~any(named('-clc')), clc, end
    
  otherwise % assume 'page'
    
    if any(ishandle(mode)), 
         fig_id = mode;
    else fig_id = gcf; varargin = [{mode} varargin];
    end    
    
    nA = sum(varargin{1} == '%')+1; 
    if any(named('-n')), nA = get_('-n'); end
    if nA == numel(varargin), varargin{end+1} = '-move'; end
    
    page = [ps_folder sprintf(varargin{1:nA})]; 
    figs_to_ps(fig_id, page, varargin{nA+1:end}); 
    
    pause(0.05), close(fig_id)
    
end







function [nFigs,fig_handles]=figs_to_ps(varargin)
% [nfigs,fighandles]=figs_to_ps([filename], [-flags])
%   prints all figures to output.ps
%   
%   nFigs: number of figures printed
%   fig_handles: handles to open figures
%   
%   filename specifies the output file name (default: 'output.ps'), 
%   which is created in both the local temp directory and the current
%   working directory (workaround for 2015a+)
%
%   Optional flags: -erase: starts from a blank file (instead of appending)
%                   -delete: synonym for -erase
%                   -portrait: print in portrait orientation
%                   -landscape: print in landscape orientation [default]
%                   -silent: suppress console output
%                   -clps: if used suppresses fixing figures for illustrator
%                   -move: if used does not copy file to current folder
%
%   History:
%       CDE 26-10-2015 adapted from figs_to_ps
%       PRM added destination pwd to copyfile
%       SP 29-10-2015 Added the clps functionality to it, default is doing
%       it, also added flag for moving file, not needed if you use
%       combine_PStoPDF
%       CDE now skips figures which are Visible = 'off'
%       PRM added -delete as synonym for -erase
%       PRM 25-Jan-2017 append .ps to filename if absent
%       PRM 12-Feb-2018 make 

hidx = cellfun(@(x)ishandle(x(1)),varargin); 
if any(hidx)
    all_handles = [varargin{hidx}];
    varargin(hidx) = [];
else all_handles = get(groot,'Children');
end


filename = find(~strncmp('-',varargin,1),1);

if isempty(filename), filename = 'output.ps'; 
else                  filename = varargin{filename};
    
end

if ~strncmp(filename(end-2:end),'.ps', 3)
    filename = [filename, '.ps'];
end

if any(strncmpi(varargin,'-erase',2) | strncmpi(varargin,'-delete',2)) ... 
                                    && exist([tempdir filename],'file')
    delete([tempdir filename])
end

volubleFlag = ~any(strncmpi(varargin,'-silent',2));

if any(strncmpi(varargin,'-portrait',2)), paper_orientation = 'portrait';
else                                      paper_orientation = 'landscape';
end

% index for cleaning post script, standard is on
if any(strncmpi(varargin,'-clps',2)), clean_figs = 0;
else                                  clean_figs = 1;
end

% index for moving (copying) file to current folder, standard is on
if any(strncmpi(varargin,'-move',2)), move_file = 0;
else                                  move_file = 1;
end

% index for resizing the figure to the default page size, standard is ON
if any(strncmpi(varargin,'-resize',2)), preserveSize = 1;
else                                    preserveSize = 0;
end

% index for resizing the figure to the default page size, standard is ON,
% SP 07-2017
if any(strncmpi(varargin,'-bestfit',2)), best_fit = 0;
else                                     best_fit = 1; %#ok<*SEPEX>
end
% This is a guess as to the version number where -bestfit got introduced
best_fit = best_fit & ~verLessThan('matlab','9.0.0'); 


nFigs=0;
if (isempty(all_handles)), return, end

[~,order] = sort(cat(1,all_handles.Number),'ascend');
all_handles = all_handles(order);

fig_handles=gobjects(0);

for ii=1:size(all_handles,1)
    if ~strcmp(all_handles(ii).Type,'figure'), continue, end
    if ~strcmp(all_handles(ii).Visible,'on'), continue, end
    
    fig_handles=[fig_handles; all_handles(ii)]; %#ok<AGROW>
    figure(all_handles(ii));

    orient(paper_orientation);
    set(gcf,'Renderer','painters')

    if preserveSize
        set(gcf,'PaperPositionMode','auto')
        print(all_handles(ii),[tempdir filename],'-dpsc2','-r0'); % '-dpsc2'
    elseif best_fit
        print(all_handles(ii),[tempdir filename],'-dpsc2','-append', '-bestfit'); % '-dpsc2'
    else
        print(all_handles(ii),[tempdir filename],'-dpsc2','-append'); % '-dpsc2'
    end

    if volubleFlag
        fprintf('%s: Saved figure %d to file %s \n', ...
                                mfilename, ...
                                all_handles(ii).Number, filename)
    end
end

% clean up figures linestyle etc (clps function)
if clean_figs
    fprintf('fixing figures linestyle\n')
    fixPSlinestyle([tempdir filename])
end

if move_file
    fprintf('moving files to pwd folder\n')
    copyfile([tempdir filename], pwd)
end

nFigs=size(fig_handles,1);

%% subfunction fixPSlinestyle %%
function fixPSlinestyle(ps_in)

% Error checking
narginchk(1, 2);
if ~ischar(ps_in) || (nargin == 2 && ~ischar(ps_out))
  error('Input arguments must be file names (char).');
end

% Open file and read it in
fid = fopen(ps_in, 'r');
str = fread(fid);
str = char(str');
fclose(fid);

% 0 setlinecap and 2 setlinecap need to be 1 setlinecap (also join, PRM
% 2015)
str = strrep(str, '0 setlinecap', '1 setlinecap 1 setlinejoin');
str = strrep(str, '2 setlinecap', '1 setlinecap 1 setlinejoin');

% Replace spurious calls to non-existent font PRM 2017
str = strrep(str, 'mwb_cmsy10', 'Helvetica');

% Write out to file
fid = fopen(ps_in, 'w');
fprintf(fid, '%s', str);
fclose(fid);


%% combine PS files into a single PDF
function combine_PStoPDF(folder_path, destination_name, destination_path)
% function that combines all ps files in a specific folder
% New plan: first clean up files with clps, then combine to pdf
% combining to a ps file give huge files
% !!! important note, make sure there are not extra temporary ps files in
% the temp folder while using this, better to create a new folder in the
% temp folder !!!
% SP 28-10-2015
%
% uses ! to execute code that works with ghostscript

if nargin < 2, destination_name = 'combined.pdf'; end

if nargin < 1, folder_path = tempdir;             end

% first go to the folder
start_path = pwd;
cd(folder_path);

% ahum function does not work on pc, has no ghost script
% execute ghost script to combine ps files to single pdf
if ismac
    ! /usr/local/bin/gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=combined.pdf *.ps
else
    % so for pc the * thing does not work, need two lines to call gs, first
    % one makes a list of files, second one creates the pdf
    % CDE - changed 12/11/15 to work woth other folder paths. 
    system(['for %i in (' folder_path '*.ps) do ECHO "%i" >> filename.lst']);
   
    ! "C:\Program Files\gs\gs9.27\bin\gswin64c.exe" -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dEncodeColorImages=false -dEncodeGrayImages=false -dEncodeMonoImages=false -sOutputFile="combined.pdf" @filename.lst
end 

% some old code used to figure stuf out on how to call ghostscript etc
% ! for %i in (C:\Users\New\AppData\Local\Temp\collect_ps\*.ps) do ECHO "%i" >> filename.lst
% ! "C:\Program Files\gs\gs9.18\bin\gswin64c.exe" -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile="combined.pdf" @filename.lst
% 
% ! for %i in (C:\Users\New\AppData\Local\Temp\collect_ps\*.ps) do "C:\Program Files\gs\gs9.18\bin\gswin64c.exe" -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile="combined.pdf"
% gs -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS=842 -dDEVICEHEIGHTPOINTS=595 -dEPSFitPage -dNOPAUSE -dBATCH -dSAFER  -sOutputFile=OUT.pdf *.eps; open OUT.pdf

% go back to the original path
cd(start_path);

try
  % move file to start folder, possibly changing the name to destination name
  if nargin < 1
    copyfile([folder_path 'combined.pdf']);    
  elseif nargin > 2
    copyfile(fullfile(folder_path, 'combined.pdf'), fullfile(destination_path, destination_name));
  else
    copyfile(fullfile(folder_path, 'combined.pdf'), fullfile(start_path, destination_name));
  end

catch C
    
    warning(C.getReport)
    cd(folder_path)
end

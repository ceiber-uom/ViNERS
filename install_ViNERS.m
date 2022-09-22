
function install_ViNERS
% This user interface should guide you through the process of installing
% ViNERS on your machine. This is a matlab repository with several external
% dependencies which is designed to work within a SPARC data structure.
% This UI will help you download and install gmsh, neuron, eidors
% (each a dependency of eidors) and should work for PC, mac, or linux. 
% 
% Tested on several different versions of matlab (2019b - 2022a), should
% work without issues if your version of matlab supports uifigure
% 
%  https://github.com/ceiber-uom/ViNERS/wiki/Installation
%  https://www.biorxiv.org/content/10.1101/2021.02.10.430563v2
% 
% Version History: 
% v1.0 - wrote it. 
% v1.3 - change to github path (formerly the UoM public gitlab)
% 
% 22 September 2022 - Calvin Eiber <c.eiber@ieee.org>

set(0,'ShowHiddenHandles','on'), delete(get(0,'Children'))
set(0,'ShowHiddenHandles','off')

splashScreen('reset'); % reset

diary('ViNERS-installation.log')
fprintf('ViNERS installation %s\n', datestr(now))

done = onCleanup(@() diary('off'));

if ismac, setup_ViNERS_on_mac
elseif isunix, setup_ViNERS_on_unix    
elseif ispc, setup_ViNERS_on_PC
else error('Platform not supported, please install manually')
end


function setup_ViNERS_on_unix

h = waitbar(0, 'Gathering required ion ...');
stop_if_cancelled = @() assert(ishandle(h), 'cancelled by user'); 
p_ = @(x, varargin) [x.folder filesep x.name varargin{:}];
s_ = @(varargin) system(sprintf(varargin{:}));

n_steps = 5;
[not_found,out] = system('which nrniv'); % must be on path

config.neuron = firstItem(out,not_found);
stop_if_cancelled(); h = waitbar(1/n_steps, h);

%% EIDORS

out = which('eidors_startup.m');
if isempty(out)
    
  if isempty(userpath)
    if exist('~/MATLAB/','dir')
       out = dir('~/MATLAB/**/eidors_startup.m'); 
    end
  else out = dir([userpath '/**/eidors_startup.m']);
  end
  if ~isempty(out), out = (out(1)); end
end
config.eidors = firstItem(out,isempty(out));
if ~isempty(config.eidors)
  config.eidors = strrep(config.eidors,'eidors_startup.m','');
end

stop_if_cancelled(); h = waitbar(2/n_steps, h);

%% GMSH

[not_found,out] = system('which gmsh');
if not_found % check pip
   [not_found,out] = system('pip show gmsh');
   if not_found
    out = regexp(out,'(?<=Location: )[^\n]*','match','once');
    out = dir([out '/**/gmsh']);
    out = out(strcmp({out.name},'gmsh'));
    if ~isempty(out), out = p_(out(1)); not_found = 0; end
   end
end

config.gmsh = firstItem(out,not_found);
stop_if_cancelled(); h = waitbar(3/n_steps, h);

%% GIT

[not_found,~] = system('which git');
if not_found, config.git = [];    
else config.git = 'git';
end

config.install_path = [pwd filesep 'ViNERS']; 

stop_if_cancelled(); h = waitbar(4/n_steps, h);

%% Parallel
config.needs_PCT  = isempty(which('getCurrentTask'));

%% LOGO (also helps determine what kind of setup to do)

if exist('./source-data/schematic.png','file')
 config.img = imread('./source-data/schematic.png','backgroundColor',[1 1 1]);
 config.img = double(config.img)/255;
 config.state = 'welcome';
 config.install_type = 'full';
elseif exist('../source/schematic.png','file')   
 config.img = imread('../source/schematic.png','backgroundColor',[1 1 1]);
 config.img = double(config.img)/255;
 config.state = 'welcome-repair';
 config.install_type = 'partial';
else
 config.img = ones(10,10,3);
 config.state = 'welcome-download';
 config.install_type = 'full';
end

stop_if_cancelled(); h = waitbar(5/n_steps, h);
pause(0.2), delete(h);

%% Begin interactive configuration

while ~strcmp(config.state,'done')

  old_state = config.state;
    
  [config,handles] = splashScreen(config);
  
  switch old_state
    case {'welcome','welcome-download'}
      config.install_path = handles(6).Value;
      if ~exist(config.install_path,'dir'), mkdir(config.install_path); end
      cd(config.install_path), pause(0.05)
    case 'get-NRN',    config.neuron = handles(6).Value;
    case 'get-EIDORS', config.eidors = handles(6).Value;
    case 'get-GMSH',   config.gmsh   = handles(6).Value;
  end
    
  switch config.state
   case 'git pull viners' % download 
       
    url = gitlab_path;
    branch = 'ViNERS-master';
    
    if exist(config.install_path,'dir'), cd '..', pause(0.05)        
      list = dir(config.install_path);
      list(cellfun(@(n) all(n=='.'),{list.name})) = [];
      if ~isempty(list)
        sel = uiconfirm(handles(1).Parent, sprintf( ... 
                  'Erase %s (%d files, %d folders)?', config.install_path, ...
                  sum(~[list.isdir]),sum([list.isdir])), 'Confirm download path');
        if strcmp(sel,'Cancel'), config.state = old_state; continue, end
      end
         rmdir(config.install_path,'s');
    end, mkdir(config.install_path); cd(config.install_path); pause(0.05)
    
    if ~isempty(config.git)
        s_('%s clone --depth=1 %s "%s/%s"',config.git,url,config.install_path,branch);
    else
     url = [url '/-/archive/master/ViNERS-master.tar.gz']; %#ok<AGROW>

     % web(url,'-browser')
     try 

      u = webread(url); 
      fi = fopen([tempdir '/install-Viners.tar.gz'],'w');
      fwrite(fi,u,class(u));
      fclose(fi);
      untar([config.install_path filesep branch '.tar.gz'],config.install_path);
      % untar([tempdir '/install-Viners.tar.gz'],config.install_path);
      delete([config.install_path '/pax_global_header'])
     catch C %#ok<NASGU>
      url = strrep(url,'.tar.gz','.zip');
      web(url,'-browser')
      config = splashScreen(config,'manual-unzip');

      list = dir([config.install_path '/*.zip']);
      if isempty(list)
          error('Please download ViNERS from <a href="%s">GitLab</a> and run the contained installer', url)
      end

      list = list([list.datenum] == max([list.datenum]));
      s_('unzip -f -o "%s"', p_(list));
     end
    end
    
    
    s_('mv -f %s/%s/{.,}* %s',config.install_path,branch,config.install_path)
    rmdir([config.install_path filesep branch],'s')
    
    config.state = 'list-todo';
   case 'do-unpack'
       
       
       try tools.file('info')
       catch C % probably missing the compiled semaphore code e.g. unix
           cd './+tools'
           ! mex -O -v semaphore.c
           cd ..
       end
       
       % path('./code',path)
       tools.make_SPARC_structure      
       
       if exist('./code','dir'), cd code
       elseif exist('../code','dir'), cd ../code % probably does nothing
       else warning('Could not find CODE folder. ViNERS may not be configured correctly.')
       end
       
       % tools.configuration('make-new-file',config); % do after nrnivmodl
       config.state = 'notify-nrnivmodl';
        
       if exist('./source-data','dir')
           system('mv source-data source');
           rmdir('../source','s') % remove empty folder source folder
           system('mv -f source ..'); 
       else 
         warning('Could not find "source-data" folder. ViNERS may not be configured correctly.')
       end
       
   case 'do-nrnivmodl'       
      
       cd ../source/nrnmech/
       
       clc
     ! sudo apt install -y make libopenmpi-dev libmeschach-dev       
     ! nrnivmodl
       clc
       
       u = dir('./*/special');
       u = u([u.datenum] == max([u.datenum]));
       config.neuron = p_(u(1));
       
       cd ../../code
       
       tools.configuration('make-new-file',config);
       config.state = 'congrats';

   case {'list-todo','get-NRN','get-EIDORS','get-GMSH','notify-unpack'}, continue % no processing to do
       
   otherwise
     disp(['unknown state "' config.state '", ignoring.'])
  end
end
%%

fprintf('\n%s\n\n      ViNERS installation complete\n\n%s\n\n', ...
                '='*ones(1,40), '='*ones(1,40))
delete(handles(1).Parent)            

function setup_ViNERS_on_mac
error('Setting up ViNERS on MAC is left as an exerise to the mac user.')

function setup_ViNERS_on_PC

h = waitbar(0, 'Gathering required information ...');
stop_if_cancelled = @() assert(ishandle(h), 'cancelled by user'); 
p_ = @(x, varargin) [x.folder filesep x.name varargin{:}];
s_ = @(varargin) system(sprintf(varargin{:}));

n_steps = 5;

[~,pf] = system('echo %programfiles%');
pf = [strtrim(pf) filesep];

[not_found,out] = system('where nrniv.exe'); % must be on path
if not_found, out = dir('C:/n*/**/nrniv.exe');
  if ~isempty(out), out = p_(out(1)); not_found = 0; end
end

config.neuron = firstItem(out,not_found);
stop_if_cancelled(); h = waitbar(1/n_steps, h);

%% EIDORS

out = which('eidors_startup.m');
if isempty(out)
  out = dir([userpath '/**/eidors_startup.m']); 
  if ~isempty(out), out = (out(1)); end
end
config.eidors = firstItem(out,isempty(out));
if ~isempty(config.eidors)
  config.eidors = strrep(config.eidors,'eidors_startup.m','');
end

stop_if_cancelled(); h = waitbar(2/n_steps, h);

%% GMSH

[not_found,out] = system('where gmsh.exe');
if not_found, out = dir([pf '/G*/**/gmsh.exe']); 
  if ~isempty(out), out = p_(out(1)); not_found = 0; end
end

config.gmsh = firstItem(out,not_found);
stop_if_cancelled(); h = waitbar(3/n_steps, h);

%% GIT

[not_found,~] = system('where git');
if not_found, config.git = [];    
else config.git = 'git';
end

config.install_path = [pwd filesep 'ViNERS']; 
stop_if_cancelled(); h = waitbar(4/n_steps, h);

%% Parallel
config.needs_PCT  = isempty(which('getCurrentTask'));

%% LOGO (also helps determine what kind of setup to do)

if exist('./source-data/schematic.png','file')
 config.img = imread('./source-data/schematic.png','backgroundColor',[1 1 1]);
 config.img = double(config.img)/255;
 config.state = 'welcome';
 config.install_type = 'full';
elseif exist('../source/schematic.png','file')   
 config.img = imread('../source/schematic.png','backgroundColor',[1 1 1]);
 config.img = double(config.img)/255;
 config.state = 'welcome-repair';
 config.install_type = 'partial';
else
 config.img = ones(10,10,3);
 config.state = 'welcome-download';
 config.install_type = 'full';
end

stop_if_cancelled(); h = waitbar(5/n_steps, h);
pause(0.05); 

%% Begin interactive configuration

path(path,pwd)

add_ = @(a,b) ~strncmpi(fliplr(a),fliplr(b),length(b)) && ...
                      exist([a filesep b],'file');

while ~strcmp(config.state,'done')

  old_state = config.state;
    
  [config,handles] = splashScreen(config);
  
  if ishandle(h), delete(h); end
  
  switch old_state
    case {'welcome','welcome-download'}
      config.install_path = handles(6).Value;
      
      if any(config.install_path == ' ')
      
        sel = uiconfirm(handles(1).Parent,['NEURON may not be able to ' ... 
                   'correctly handle file-paths with spaces, continue?'], ...
                   'Confirm install path');
        if strcmp(sel,'Cancel'), config.state = old_state; continue,    end
      end   
      if ~exist(config.install_path,'dir'), mkdir(config.install_path), end
      cd(config.install_path), pause(0.05)
    case 'get-NRN',    config.neuron = handles(6).Value;        
        if add_(config.neuron,'\nrniv.exe')
            config.neuron = fullfile(config.neuron,'nrniv.exe'); 
        end
    case 'get-EIDORS', config.eidors = handles(6).Value;        
    case 'get-GMSH',   config.gmsh   = handles(6).Value;
      if add_(config.gmsh,'\gmsh.exe')
          config.gmsh = fullfile(config.gmsh,'gmsh.exe'); 
      end
  end
    
  switch config.state
   case 'git pull viners' % download 
       

    url = gitlab_path; 
    branch = 'ViNERS-master';
    
    if exist(config.install_path,'dir'), cd '..', pause(0.05)        
      list = dir(config.install_path);
      list(cellfun(@(n) all(n=='.'),{list.name})) = [];
      if ~isempty(list)
        sel = uiconfirm(handles(1).Parent, sprintf( ... 
                  'Erase %s (%d files, %d folders)?', config.install_path, ...
                  sum(~[list.isdir]),sum([list.isdir])), 'Confirm download path');
        if strcmp(sel,'Cancel'), config.state = old_state; continue, end
      end
         rmdir(config.install_path,'s');
    end, mkdir(config.install_path); cd(config.install_path); pause(0.05)
    
    if ~isempty(config.git)
        s_('%s clone --depth=1 %s "%s/%s"',config.git,url,config.install_path,branch);
    else
     url = [url '/-/archive/master/ViNERS-master.tar.gz']; %#ok<AGROW>

     % web(url,'-browser')
     try 

      u = webread(url); 
      fi = fopen([tempdir '/install-Viners.tar.gz'],'w');
      fwrite(fi,u,class(u));
      fclose(fi);
      untar([config.install_path filesep branch '.tar.gz'],config.install_path);
      % untar([tempdir '/install-Viners.tar.gz'],config.install_path);
      delete([config.install_path '/pax_global_header'])
     catch C %#ok<NASGU>
      url = strrep(url,'.tar.gz','.zip');
      web(url,'-browser')
      config = splashScreen(config,'manual-unzip');

      list = dir([config.install_path '/*.zip']);
      if isempty(list)
          error('Please download ViNERS from <a href="%s">GitLab</a> and run the contained installer', url)
      end

      list = list([list.datenum] == max([list.datenum]));
      s_('unzip -f -o "%s"', p_(list));
     end
    end
    
    s_('for /d %%A in ("%s\\%s\\*") do move "%%A" "%%~dpA.."',config.install_path,branch);
    s_('MOVE "%s\\%s\\*" "%s"',config.install_path,branch,config.install_path);
    rmdir([config.install_path filesep branch],'s')
    
    config.state = 'list-todo';
   case 'do-unpack'
       
       % path('./code',path)
       tools.make_SPARC_structure      
       
       if exist('./code','dir'), cd code
       elseif exist('../code','dir'), cd ../code % probably does nothing
       else warning('Could not find CODE folder. ViNERS may not be configured correctly.')
       end
       
       tools.configuration('make-new-file',config);
       config.state = 'manual-mknrndll';
        
       if exist('./source-data','dir')           
         system('ren source-data source');
         rmdir('../source','s') % remove empty folder source folder
         system('move /Y source ..'); 
       else 
         warning('Could not find "source-data" folder. ViNERS may not be configured correctly.')
       end
       
   case 'move-nrndll'       

     obj = dir('../source/nrnmech/*.dll');       
     s_('copy "%s" ".\\%s"',p_(obj(1)),obj(1).name);   
     config.state = 'congrats';

   case {'list-todo','get-NRN','get-EIDORS','get-GMSH','notify-unpack','done'}, continue % no processing to do
       
   otherwise
     disp(['unknown state "' config.state '", ignoring.'])
  end
end
%%

fprintf('\n%s\n\n      ViNERS installation complete\n\n%s\n\n', ...
                '='*ones(1,40), '='*ones(1,40))
delete(handles(1).Parent)

%%
 
% // example CONFIGURATION file generated 11-Aug-2020 02:30:55
% 
% { "env-list": [{
%     "name": "THE-STREAM",
%     "user": "Calvin",
%     "machine": "x64-based PC",
%     "root": "C:\Users\Calvin\Documents\MATLAB\Keast-lab\oSPARC"
%     "cache": "pn-mdl-%d/",
%     "octave": "false", 
%     "gmsh": "C:\Program Files\GMSH 4.8.3\gmsh.exe",
%     "neuron": "C:\nrn\bin\nrniv.exe",
%     "eidors": "C:\Users\Calvin\Documents\MATLAB\Keast-lab\EIDORS\eidors",
%     "toolbox-gsa": "null"
%   }]
% }


function [config,handles] = splashScreen(config,command,varargin)

persistent h
if ischar(config), h = []; return, end
if isempty(h)
    %%
    h = uifigure('Position',[300 340 780 600]);
    uiaxes('Parent',h,'Position',[400 60 360 480],'Color','none');
    imagesc((config.img + 0.3)/1.3,'Parent',h.Children(end));
    axis(h.Children,'image','off');
    
    t = repmat('instructionText',[1 1000]);
    t(1:100:end) = newline;
    t(1:20) = 'A';
    
    %%
    delete(h.Children(1:end-1))
    uilabel(h,'Text','ViNERS setup','Position',[15 520 360 20],'FontSize',16,'VerticalAlignment','top','FontWeight','bold');
    uilabel(h,'Text',t,'Position',[15 300 660 210],'FontSize',14,'VerticalAlignment','top');
    
    uibutton(h,'push','Text','Button','Position',[15 180 360 26],'fontSize',12);
    
    uidropdown(h,'Items',{'option 1','option2'},'Position',[15 260 360 26],'fontSize',12);
    uieditfield(h,'text','Value','textedit','Position',[15 220 360 26],'fontSize',12);
    uibutton(h,'push','Text','Button','Position',[15 180 360 26],'fontSize',12);
    uibutton(h,'push','Text','Link','Position',[15 140 360 26],'fontSize',12);
    
    uibutton(h,'push','Text','Exit','Position',[15 15 120 40],'fontSize',14,'BackgroundColor',[1 0.85 0.85]);
    uibutton(h,'push','Text','Back','Position',[150 15 120 40],'fontSize',14);    
    uibutton(h,'push','Text','Next','Position',[400 15 360 40],'fontSize',14,'BackgroundColor',[0.85 1 0.85]);
    
    h.Children(end).UserData = numel(h.Children);
   %% 
else assert(ishandle(h),'cancelled by user')
    delete(h.Children(1:end-h.Children(end).UserData))
end

h.Children(2).Visible = 'on'; % back button always visible unless...
h.Children(6).Editable = true; % text editable by default 

h.Children(8).Visible = 'off'; % sudo apt-get
h.Children(7).Visible = 'off'; % dropdown
h.Children(6).Visible = 'off'; % textedit
h.Children(5).Visible = 'off'; % act button
h.Children(4).Visible = 'off'; % link button

if nargin == 1, command = config.state; end




switch command
  case {'welcome','welcome-download'}
    %%
    h.Children(9).Text = ['Welcome to the setup for the ViNERS ' newline ...
                          '(Visceral Nerve Ensemble Recording & Stimulation)' newline ...
                          'computational pipeline. ' newline newline ...
                          'This helper will ensure that you have all of the required ' newline ...
                          'components for ViNERS installed on your computer ' newline ...
                          'and set up ViNERS for your machine.' newline newline ...
                          'to begin, please select a directory where ' newline ...
                          'ViNERS will be downloaded and installed.']; 
    
    h.Children(6).Value = config.install_path;
    h.Children(6).Position = [15 260 360 26];
    h.Children(6).Visible = 'on';  % textedit
    
    h.Children(5).Text = 'browse';
    h.Children(5).Position = [285 225 90 26];
    h.Children(5).ButtonPushedFcn = @(~,~) browseFile(h.Children(6),'Select installation directory');
    h.Children(5).Visible = 'on';  % button
        
    state = {'list-todo','welcome'};
    
    if strcmp(command,'welcome-download'), state{1} = 'git pull viners'; end
    % return
  case 'welcome-repair'
      
    h.Children(9).Text = ['Welcome to the setup for the ViNERS ' newline ...
                  '(Visceral Nerve Ensemble Recording & Stimulation)' newline ...
                  'computational pipeline. ' newline newline ...
                  'Please select from one of the following options:']; 

    state = {'get-NRN','get-EIDORS','get-GMSH','manual-mknrndll','do-unpack','welcome-download'};
    
    h.Children(7).Items = {'reinstall NEURON','reinstall EIDORS','reinstall GMSH', ...
                           'generate nrnmech.dll','run tools.make_SPARC_structure', ...
                           'download and reinstall ViNERS'};
    h.Children(7).Position(2) = 390;
    h.Children(7).Visible = 'on'; 
    
    
    h.UserData = 0;

    for b = 1:3
        h.Children(b).ButtonPushedFcn = @(~,~) set(h,'UserData',b);
    end

    while h.UserData == 0, pause(0.02), end
    if h.UserData == 3, delete(h), error('cancelled'), end
    if h.UserData == 2, config.state = 'welcome-repair'; end
    if h.UserData == 1
        idx = strcmp(h.Children(7).Items,h.Children(7).Value);
        config.state = state{idx}; 
    end
    
    handles = h.Children; 
    h.Children(3).ButtonPushedFcn = @(~,~) close(h); % exit button

    return
      
  case 'manual-unzip'
            
    h.Children(9).Text = ['Automatic download failed, please download ViNERS from ' newline ...
                          'GitLab into ' config.install_path newline ...
                          'I recommend cancelling this installer and running ' newline ...
                          'the installer packaged in the downloaded .zip file.'];
    h.Children(8).Position(3) = 600;
    
    h.Children(5).Text = 'open ViNERS on GitLab';
    h.Children(5).Position = [15 425 360 26];
    h.Children(5).ButtonPushedFcn = @(~,~) web(gitlab_path,'-browser');
    h.Children(5).Visible = 'on';  % button
    
    
    state = {'manual-unzip','welcome'};
    
  case 'list-todo'
              
      lost = ''; 
      found = ''; 
      
      if isempty(config.neuron), lost = [lost newline 'NEURON (neuron.yale.edu)'];
      else found = [found newline 'NEURON (neuron.yale.edu)'];
      end
      
      if isempty(config.eidors), lost = [lost newline 'EIDORS (eidors.org)'];
      else found = [found newline 'EIDORS (eidors.org)'];
      end
      
      if isempty(config.gmsh), lost = [lost newline 'GMSH (gmsh.info)'];
      else found = [found newline 'GMSH (gmsh.info)'];
      end
      
      if config.needs_PCT, lost = [lost newline 'MATLAB Parallel Computing Toolbox']; end
      
      if ~isempty(lost), lost = [newline newline 'The following were not located on your system. ' newline ...
                                   'This may be because they have yet to be installed, ' newline ...
                                   'or were installed to an unusual location:' lost newline]; end
      if ~isempty(found), found = ['The following were already located on your system:' found]; end
      
      h.Children(9).Text = [found lost];

      state = {'get-NRN','welcome'};
      
    case 'get-NRN'
      
      state = {'get-EIDORS','list-todo'};
      if isempty(config.neuron), state{1} = 'get-NRN';
        h.Children(9).Text = ['No installation of NEURON was found on your' newline ...
                              'system. If you have already installed NEURON, please ' newline ...
                              'enter the path to your NEURON installation. If you ' newline ... 
                              'have not yet installed neuron, you should download ' newline ...
                              'NEURON using the link below and install it now.' newline newline ...
                              'When you have done so, please enter the path to ' newline ...
                              'NEURON (nrniv.exe):'];
        config.neuron = '/path/to/neuron';
      else
        h.Children(9).Text = ['An installation of NEURON was found on your' newline ...
                              'system at ' config.neuron newline newline ... 
                              'If this is OK, press continue. Otherwise, enter the path ' newline ...
                              'to the NEURON installation you wish to use with ViNERS' newline ...
                              '(nrniv.exe):'];
      end
      
      url = 'https://neuron.yale.edu/neuron/download';
      if ispc, nrn_tgt = 'nrniv.exe'; 
      else     nrn_tgt = 'nrniv';
      end
      
      h.Children(7).Visible = 'off'; % dropdown
      
      h.Children(6).Value = config.neuron;
      h.Children(6).Position = [15 350 360 26];
      h.Children(6).Visible = 'on';  % textedit

      h.Children(5).Text = 'browse';
      h.Children(5).Position = [285 318 90 26];
      h.Children(5).ButtonPushedFcn = @(~,~) browseFile(h.Children(6),'Select NEURON directory',nrn_tgt);
      h.Children(5).Visible = 'on';  % button

      h.Children(4).Text = 'NEURON website (installer)';
      h.Children(4).Position = [15 318 262 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button      
      
      if isunix
          h.Children(8).Text = 'sudo apt install neuron';
          h.Children(8).Position = [15 320 262 26];
          h.Children(8).ButtonPushedFcn = @(b,~) apt_install(b,'neuron neuron-dev');
          h.Children(8).Visible = 'on';  % button      
      end
      
    case 'get-EIDORS'
        
     
      state = {'get-GMSH','get-NRN'};
      
      if isempty(config.eidors), state{1} = 'get-EIDORS';
        h.Children(9).Text = ['The EIDORS toolbox for MATLAB was not found on your' newline ...
                              'system. If you have already downloaded EIDORS, please ' newline ...
                              'enter the path to your EIDORS toolbox (specifically, the ' newline ...
                              '"eidors_startup.m" file, located in the /eidors/ directory).' newline ... 
                              'If you have not yet downloaded EIDORS, you should download ' newline ...
                              'EIDORS using the link below.' newline newline ...
                              'When you have done so, please enter the path to' newline ...
                              'EIDORS (eidors_startup.m):'];
        config.eidors = '/path/to/eidors';
      else
        h.Children(9).Text = ['The EIDORS toolbox was found on your' newline ...
                              'system at ' config.eidors newline newline ... 
                              'If this is OK, press continue. Otherwise, enter the path ' newline ...
                              'to the EIDORS installation you wish to use with ViNERS' newline ...
                              '(eidors_startup.m):'];
      end
      
      url = 'https://sourceforge.net/projects/eidors3d/';
      eidors_tgt = 'eidors_startup.m';
      
      h.Children(6).Value = config.eidors;
      h.Children(6).Position = [15 350 360 26];
      h.Children(6).Visible = 'on';  % textedit

      h.Children(5).Text = 'browse';
      h.Children(5).Position = [285 318 90 26];
      h.Children(5).ButtonPushedFcn = @(~,~) browseFile(h.Children(6),'Select EIDORS directory',eidors_tgt);
      h.Children(5).Visible = 'on';  % button

      h.Children(4).Text = 'EIDORS toolbox (download)';
      h.Children(4).Position = [15 318 262 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button      
      
    case 'get-GMSH'
             
      state = {'notify-unpack','get-EIDORS'};
      
      if config.needs_PCT, state{1} = 'get-PCT'; end
      
      if isempty(config.gmsh), state{1} = 'get-GMSH';
        h.Children(9).Text = ['No installation of GMSH was found on your' newline ...
                              'system. If you have already installed GMSH, please ' newline ...
                              'enter the path to your GMSH installation. If you ' newline ... 
                              'have not yet installed neuron, you should download ' newline ...
                              'GMSH using the link below and install it now.' newline newline ...
                              'When you have done so, please enter the path to ' newline ...
                              'GMSH (gmsh.exe):'];
        config.gmsh = '/path/to/gmsh';
      else
        h.Children(9).Text = ['An installation of GMSH was found on your' newline ...
                              'system at ' config.gmsh newline newline ... 
                              'If this is OK, press continue. Otherwise, enter the path ' newline ...
                              'to the GMSH installation you wish to use with ViNERS' newline ...
                              '(gmsh.exe):'];
      end
      
      url = 'https://gmsh.info';
      
      if ispc, gmsh_exe = 'gmsh.exe'; 
      else     gmsh_exe = 'gmsh';
      end
      
      
      h.Children(6).Value = config.gmsh;
      h.Children(6).Position = [15 350 360 26];
      h.Children(6).Visible = 'on';  % textedit

      h.Children(5).Text = 'browse';
      h.Children(5).Position = [285 318 90 26];
      h.Children(5).ButtonPushedFcn = @(~,~) browseFile(h.Children(6),'Select GMSH',gmsh_exe);      
      h.Children(5).Visible = 'on';  % button

      h.Children(4).Text = 'GMSH website (binary)';
      h.Children(4).Position = [15 318 262 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');
      h.Children(4).Visible = 'on';  % button
      
      
      if isunix
          h.Children(8).Text = 'pip install gmsh';
          h.Children(8).Position = [15 320 262 26];
          h.Children(8).ButtonPushedFcn = @(b,~) apt_install(b,'gmsh');
          h.Children(8).Visible = 'on';  % button      
      end
      
    case 'get-PCT'
        
        state = {'notify-unpack','get-GMSH'};

        h.Children(9).Text = ['ViNERS uses the MATLAB Parallel Computing Toolbox' newline ...
                              'to accelerate many of its computations. ' newline newline ...
                              'While strictly speaking you will be able to use ViNERS' newline ... 
                              'without the Parallel Computing Toolbox, you may need to' newline ...
                              'modify some functions and large simulations will be slow.' newline newline ...
                              'After installation is complete, in the toolstrip, please go to:' newline ...
                              'Home tab -> Environment -> Add-Ons - > get Add-Ons, ' newline ...
                              'then search for the Parallel Computing Toolbox,' newline ...
                              'then click install. This will restart MATLAB.' ];
                          
      url = 'https://au.mathworks.com/products/parallel-computing.html';
      
      h.Children(4).Text = 'Parallel Computing Toolbox (MathWorks.com)';
      h.Children(4).Position = [15 355 262 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button
      
      
    case 'notify-unpack'
             
      state = {'do-unpack','get-GMSH'};
            
      h.Children(9).Text = ['The dependencies of ViNERS have now been installed:' newline ...
                            config.neuron newline config.eidors newline config.gmsh newline newline ...      
                            'ViNERS was designed to work within a SPARC' newline ... 
                            'data structure (roughly based on the BIDS standard). ' newline ...
                            'ViNERS will now create /primary, /code, and /source in' newline ...
                            'the installation directory:' ];
      
      h.Children(6).Visible = 'on'; % textedit
      h.Children(6).Position = [15 330 360 26];
      h.Children(6).Value = config.install_path; 
      h.Children(6).Editable = false; % static text but copyable


      url = 'https://sparc.science/help/3FXikFXC8shPRd8xZqhjVT';

      h.Children(4).Text = 'more information about the SPARC data structure';
      h.Children(4).Position = [15 255 360 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button
      
    case 'notify-nrnivmodl'
        
      state = {'do-nrnivmodl','cant-undo-unpack'};
            
      h.Children(9).Text = ['The final step is to compile the .mod files needed' newline ...
                            'by NEURON to simulate non-standard ion channels.' newline newline ...                            
                            'On linux, this can be done automatically by nrnivmodl.' newline ...
                            'You may need to enter information in at the MATLAB ' newline ...
                            'console if compiler components require installation.'];
      url = 'https://neuron.yale.edu/neuron';
      
      h.Children(4).Text = 'more information about NEURON and nrnivmodl';
      h.Children(4).Position = [15 320 360 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button


    case 'manual-mknrndll'
             
      state = {'move-nrndll','cant-undo-unpack'};
            
      h.Children(9).Text = ['The final step is to compile the .mod files needed' newline ...
                            'by NEURON to simulate non-standard ion channels.' newline ...
                            'Unfortunately, I couldn''t automate this step completely.' newline ...
                            'However, I can provide step-by-step instructions:' newline newline ...
                            '1: open a bash shell (bash installs with NEURON)' newline ...
                            '     by typing ''bash'' in the start menu.' newline newline ...                            
                            '2: type (or copy and paste by right-clicking) and press enter: ' newline newline newline newline ...
                            '3: type "mknrndll" and press enter to compile ' newline ...
                            '     the .mod files packaged with ViNERS. '];
      h.Children(8).Position([2 4]) = [200 310];
      url = 'https://neuron.yale.edu/neuron';

      dll_path = [config.install_path '/source/nrnmech']; % where I think this is
      if ~exist(dll_path,'dir')
          dll_path = strrep(dll_path','code/source','source'); % reinstall
      end
            
      h.Children(6).Visible = 'on'; % textedit
      h.Children(6).Position = [15 330 360 26];
      h.Children(6).Value = ['cd ' strrep(strrep(dll_path,'\','/'),' ','\ ')]; % formatting for bash
      h.Children(6).Editable = false; % static text but copyable
      
      h.Children(4).Text = 'more information about NEURON and mknrndll';
      h.Children(4).Position = [15 255 360 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button
      
      % Block NEXT button 
      h.Children(1).Text = 'Waiting for mknrndll...';
      h.Children(1).Enable = false;       
      dll_list = dir([h.Children(6).Value(4:end) filesep '*.dll']);
      while isempty(dll_list), pause(0.2),
            dll_list = dir([h.Children(6).Value(4:end) filesep '*.dll']);
      end
      h.Children(1).Text = 'Next';
      h.Children(1).Enable = true; 

      
    case 'cant-undo-unpack'
      
      state = {'manual-mknrndll','welcome-download'};
        
      h.Children(9).Text = ['The dependencies of ViNERS have now been installed:' newline ...
                            config.neuron newline config.eidors newline config.gmsh newline newline ...
                            'I don''t have a good way of undoing that last' newline ...
                            'step. Sorry. If you desperately need to undo,' newline ...
                            'may I recommend restarting?'];
      if isunix, state{1} = 'notify-nrnivmodl'; end
      
    case 'congrats'
      
      state = {'done','manual-mknrndll'};
      if isunix, state{2} = 'notify-nrnivmodl'; end

      h.Children(9).Text = ['Congratulations! You have now set up ViNERS on ' newline ...
                            'your computer. To get started, check out the wiki ' newline ...
                            'and have a look at the examples. Good luck and ' newline ...
                            'may your simulations reflect reality.'];
          
      url = gitlab_path('/-/wikis/home');

      h.Children(4).Text = 'more information about ViNERS';
      h.Children(4).Position = [15 255 360 26];
      h.Children(4).ButtonPushedFcn = @(~,~) web(url,'-browser');      
      h.Children(4).Visible = 'on';  % button
      
    otherwise 
      error('unknown state "%s"', config.state)
end



h.UserData = 0;

for b = 1:3
    h.Children(b).ButtonPushedFcn = @(~,~) set(h,'UserData',b);
end

while h.UserData == 0, pause(0.02), end
if h.UserData == 3, delete(h), error('Cancelled by User'), end

config.state = state{h.UserData}; 
handles = h.Children; 

h.Children(3).ButtonPushedFcn = @(~,~) close(h); % exit button


function s = firstItem(s,fail)

if nargin>1 && fail, s = []; return, end
if isstruct(s), s = [s.folder filesep s.name];
else s = strsplit(s,newline); s = s{1};
end

function browseFile(h,label,stub)
if nargin == 2, x = uigetdir(h.Value,label);
else [~,x] = uigetfile(['*' stub],label,h.Value); 
end
if ischar(x), h.Value = x; end

function apt_install(h,command)

t = h.Text; 
h.Text = 'please complete install in the matlab console';


if strcmp(command,'gmsh')
    
    system('sudo apt install python-pip');
    r = system('sudo pip install --upgrade gmsh');
    
    [~,s] = system('pip show gmsh');
    s = regexp(s,'(?<=Location: )[^\n]*','match','once');
    g = dir([s '/**/gmsh']);
    g = g(strcmp({g.name},'gmsh'));
    
    if ~isempty(g)
      h.Parent.Children(6).Value = [g(1).folder filesep g(1).name];
    end

else

    r = system(['sudo apt install ' command]);

end
if r == 0, clc, disp(['installation of ' command ' complete.']), end

h.Text = t; 

if strcmp(command,'neuron')
     [not_found,out] = system('which nrniv'); % must be on path
else [not_found,out] = system(['which ' command]); % must be on path
end

out = firstItem(out,not_found);
if ~isempty(out)
    h.Parent.Children(6).Value = out;
end

return

function p = gitlab_path(varargin)
% p = 'https://gitlab.unimelb.edu.au/lab-keast-osborne-release/ViNERS';
p = 'https://github.com/ceiber-uom/ViNERS';
p = strcat(p,varargin{:});

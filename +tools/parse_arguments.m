


function varargout = parse_arguments(user_args,varargin)
% [outputs] = parse_arguments( user_args, request_args )
% 
% Parse user varargin to get input files for model. This code, or variants 
%   of this code, was being replicated near the top of nearly every +model
%   and +plot function and could be efficiently factored out. This also
%   results in a more consistent set of command flags across functions. 
% 
% NOTE that as request_args are part of the ViNERS software, these are
% case-sensitive. 
% 
% Example: 
% [EM,AX] = tools.parse_arguments(varargin, 'LOAD', ...
%                          'eidors','eidors~/sens*.mat', 'axons');
% 
% The most efficient way to see how this is used is to copy the use from a
% working location where it is used to extract data of interest to you. 
% 
% Implemented request modes:
% ------------------------------------------------------------------------
% 'eidors' [filespec] 
%   +models accept eidors_file as the first positional argument, 
%                              as -eidors [file], or 
%                              as -file [file] (for compatibility)
%   if missing, find using uigetfile(filespec). 
% ------------------------------------------------------------------------
% 'axon' 
%   +models accept axons_file as the second positional argument, or 
%                             as -axons [file]. 
% 
%   This request also returns the membrane currents folder structure, 
%   which can be specified (in descending priority order) as any of the
%   following flags: -current, -membrane, -root, -axon-folder,
%   -axons-folder (these last three are not recommended as they may clash
%   with -axon __). the axons folder can also be inferred from the axon
%   filename e.g. ~/axons/axons (pelvic-nerve).mat defaults to:
%                 ~/axons/pelvic-nerve/ ... 
% 
% The following additonal options are specified: 
%   -anat [data] : override the .nerve field of the axon anatomy with the
%                  .nerve field in the supplied structure. 
%   -delta-xy [x y] : Displace axons by [x y] mm
% 
% ------------------------------------------------------------------------
% 'stim' 
%   +models accept stimulus_folder as the first positional argument or 
%                                  as -stim [folder].
%   if missing, use the newest folder of computed stimuli.
%   if 'eidors' also requested, the eidors filepaths will be inferred from 
%     embedded metadata if not set explicity with -eidors []. This is the
%     way in which 'stim' is actually used in practice
% 
% ------------------------------------------------------------------------
% 'wave' 
%   +models accepts waves_folder as the first positional argument or 
%                                as -waves [waves_folder].
%    (unless explicitly excluded with request arg 'w-positional-only'). 
%    waves_folder may be a full path, a tools.file input, or a name of a
%    folder in ~\waves. If missing, find using uigetdir. 
% 
% ------------------------------------------------------------------------
% 'wave+context' : Get also the EIDORS context for the waves. 
%  This is used by plots.preview_ECAP only, uses tools.INPUT_file.
% 
% ------------------------------------------------------------------------
% 'root' : In conjunction with 'axons', return also the membrane current
%          folder (usually not necessary, embedded in AX)
% 
% ------------------------------------------------------------------------
% 'fascicle'  just the fascicle anatomy requested, not finalised. 
% 
% v0.2 CDE - updated 26 Nov 2021

named = @(v) strncmpi(v,user_args,length(v)); 
get_ = @(v) user_args{find(named(v))+1};

request = @(v) strncmp(v,varargin,length(v));  % <<< CASE-SENSITIVE
par_ = @(v) varargin{find(request(v),1)+1};

ppart_ = @(s) regexp(s,'(?<=\()[^\)]+','match','once'); 
name_ = @(s) regexp(s,'(?<=[\\/])[^\\/]+$','match','once');

varargout = {}; 
opts.do_multi = any(request('multiselect')) || any(request('-m'));
opts.do_LOAD = any(request('LOAD')) || any(request('-l'));
opts.verbose = ~any(request('-q')) && ~any(request('quiet')) && ~any(named('-q'));

if opts.verbose,
     announce_load_ = @(s) fprintf('Loading %s ... \n', tools.file('T',s));
else announce_load_ = @(s) []; 
end

p_ = @(x) [x.folder filesep x.name];

if opts.do_multi, uigetfile_options = {'multiselect','on'};
else              uigetfile_options = {}; 
end

%% get EIDORS file (either stim or sens)
if any(request('eidors'))
  %%
  eidors_file = tools.file('eidors~\');
  
  % model_functions accept eidors_file as a first argument
  if numel(user_args) > 0 && ischar(user_args{1})
      eidors_file = user_args{1}; 
  end

  if any(request('stim')) && isfolder(eidors_file)    
    if any(named('-stim'))
      user_args(find(named('-stim'))+1) = {eidors_file}; 
    else
      user_args = [user_args {'-stim', eidors_file}];
    end
    named = @(v) strncmpi(v,user_args,length(v)); 
    get_ = @(v) user_args{find(named(v))+1};    
    
    % Detect eidors file from stimulus folder data
    if ~any(named('-eidors')) && ~any(named('-file'))
      try
        u = dir([eidors_file filesep '*.mat']);
        load(p_(u(1)),'notes'); %         
        if ~isempty(name_(notes{4})), notes{4} = name_(notes{4}); end
        
        eidors_file = ppart_(notes{4}); 
        if isempty(eidors_file), eidors_file = notes{4}; end

        u = par_('eidors'); 
        if any(u == '~'), u = tools.file(u); end
        u = dir(u); % what kind of eidors file requested?
        sel = contains({u.name},eidors_file); 
        if sum(sel) == 1, eidors_file = p_(u(sel)); 
        else
          n_match = sum(sel);
          sel = sel & [u.datenum] == max([u(sel).datenum]);
          warning('ViNERS:eidorsFileName',...
              '%d files in %s match "%s". Defaulting to "%s"', n_match, ...
               u(1).folder, notes{4}, u(sel).name); 
        end
      catch C, eidors_file = tools.file('eidors~\');
      end
    end
  end
  
  
  % models.function( ..., '-file', eidors_file ) also valid
  if any(named('-eidors')), eidors_file = get_('-eidors');
  elseif any(named('-file')), eidors_file = get_('-file'); 
  end
  if any(eidors_file == '~'), eidors_file = tools.file(eidors_file); end  
  
  if strncmp(eidors_file,'?',1) % user requested a menu
    
    if numel(eidors_file) == 1, eidors_file = par_('eidors'); end % '?'
    if isa(eidors_file,'function_handle'), eidors_file = eidors_file(); end    
    [eidors_file,p] = uigetfile(eidors_file(2:end), [], ...
                                 tools.file('eidors~\'), ...
                                 uigetfile_options{:});
    eidors_file = strcat(p,eidors_file); 
    
  elseif exist(eidors_file,'file') ~= 2, % use handle or menu anyway
    
    eidors_file = par_('eidors');
    if isa(eidors_file,'function_handle')
      % e.g. 'eidors',@() tools.file('get','eidors~/stim*.mat','newest')
      eidors_file = eidors_file(); % user_args{:}
    else
      if eidors_file(end) ~= filesep
        fprintf('%s not found\n', eidors_file), 
      end
     [eidors_file,p] = uigetfile(eidors_file, [], ...
                                 tools.file('eidors~\'), ...
                                 uigetfile_options{:});
      eidors_file = strcat(p,eidors_file); 
    end
  end
  
  if opts.do_multi && ~iscell(eidors_file), eidors_file = {eidors_file}; end
  if opts.do_LOAD, announce_load_(eidors_file)
    
    EM = load(eidors_file); [~,EM.filename] = fileparts(eidors_file);
    varargout = [varargout {EM}];   
  else varargout = [varargout {eidors_file}];   
  end
end

%% Find axons.mat (from models.make_population)

if any(request('axon'))
  
  axon_file = ''; 
  axons_folder = ''; 
  
  % model_functions accept axons_file/axons_folder as argument #2
  if any(named('-ax')), axon_file = get_('-ax');    
  elseif numel(user_args) > 1 && ischar(user_args{2}) && exist(user_args{2},'file')
    assert(contains(user_args{2},'axon'),'expected positional argument 2: axons')
    if isfolder(user_args{2}), axons_folder = user_args{2}; 
    else axon_file = user_args{2};
    end
  end
  
  % Get axon_folder (if not specified positionally)
  if isempty(axons_folder) || any(named('-current')) || any(named('-mem'))
    if strcmp(axon_file,'?')
      axon_file = tools.file('axons~/*.mat','prompt');
    end
    if any(named('-current')), axons_folder = get_('-current');
    elseif any(named('-mem')), axons_folder = get_('-mem');
    elseif any(named('-root')), axons_folder = get_('-root'); % provided for backwards-compatability
    elseif any(named('-axon-fo')), axons_folder = get_('-axon-fo'); 
    elseif any(named('-axons-fo')), axons_folder = get_('-axons-fo'); 
    elseif isempty(axon_file), axons_folder = ''; 
    else
      axons_folder = [tools.file('axons~/') ppart_(name_(axon_file)) '/']; 
      if ~isfolder(axons_folder), axons_folder = tools.file('axons~/'); end
    end
  end
  
  if isempty(axon_file) % determine from axons_folder 
    
    if ~isempty(axons_folder)
      axon_file = dir([axons_folder '/../a*.mat']); 
      if isempty(axon_file), axon_file = dir([axons_folder '/a*.mat']); 
      end      
    else axon_file = dir(tools.file('axons~/*.mat')); 
    end
    
    if isempty(axon_file), error('please run models.axon_poplation'), end
    label = ['(' name_(axons_folder(1:end-1)) ')'];
    sel = contains({axon_file.name}, label);
    if sum(sel) == 1, axon_file = p_(axon_file(sel));
    elseif numel(sel) == 1
      warning('ViNERS:axonFileName', ...
              '%d files in %s match %s. Defaulting to "%s"', sum(sel), ...
               axon_file(1).folder, label, axon_file(1).name); 
      axon_file = p_(axon_file(1));
    else
      n_match = sum(sel);
      if any(sel)
        sel = sel & [axon_file.datenum] == max([axon_file(sel).datenum]);
      else    sel = [axon_file.datenum] == max([axon_file.datenum]);
      end
      if sum(sel) > 1, sel = find(sel,1); end
      warning('ViNERS:axonFileName',...
              '%d files in %s match %s. Defaulting to "%s"', n_match, ...
               axon_file(1).folder, label, axon_file(sel).name); 
      axon_file = p_(axon_file(sel));
    end
  end
  if isempty(axons_folder) % derive from axons_file
    
    % get the (info) part of the filename
    ppart = regexp(axon_file,'[\\/][^\\/]+$','match','once'); % filename
    ppart = regexp(ppart,'(?<=\()[^\)]+','match','once');     % parens
    
    axons_folder = dir(regexprep(axon_file,'[\\/][^\\/]+$', ... % look here
                                          ['/' ppart '/*/index.mat'])); 
    if isempty(axons_folder), % if not found, look here next
      axons_folder = dir(regexprep(axon_file,'[\\/][^\\/]+$','/*/index.mat'));
    end
    
    if isempty(axons_folder), 
      axons_folder = tools.file('axons~/');

      warning('ViNERS:axonFileName', ...
              'Unable to find folder "%s" matching %s. Defaulting to "%s"', ...
               ppart, axon_file, axons_folder);       
    else axons_folder = [fileparts(axons_folder(1).folder) filesep]; 
    end
  end
  
  if any(request('root')), varargout = [varargout {axons_folder}]; end
  
  
  %
  if opts.do_LOAD, announce_load_(axon_file);
    
    axon_data = load(axon_file);
    
    if any(named('-anat')), axon_data.nerve = get_('-anat'); end
    if ~isfield(axon_data,'nerve') % older axons files might not have this embedded
      f_list = dir(tools.file('~/source/fascicles/*.splines.dat'));
      axon_data.nerve = mesh.read_dat_file(tools.INPUT_file(f_list, p_(axon_file)));
    end
       
    if any(named('-delta-xy'))
      dxy = get_('-delta-xy');
      fprintf('Displacing axons by [%0.3f, %0.3f] mm\n', dxy)
      axon_data.pop.fibre_xy = axon_data.pop.fibre_xy + dxy;
      axon_data.pop.unmyelinated_xy = axon_data.pop.unmyelinated_xy + dxy;
    end
    
    if ~any(named('-no-f-f')) && ~isfield(axon_data.nerve,'fascicles')
      axon_data.nerve.fascicles = axon_data.nerve.outline; 
      assert(~iscell(axon_data.nerve.fascicles));
    end
    
    axon_data.folder = axons_folder; 
    [~,axon_data.name] = fileparts(axon_file); 
    
    % stub = regexp(axon_data.name,'(?<=\()[^\)]+','match','once');
    % if ~isempty(stub) && exist([axon_data.folder stub],'dir')
    %     axon_data.name = stub; 
    %     % We now do this more transparently for everything, above. 
    %     % axon_data.folder = [axon_data.folder stub filesep];
    % end
    
       varargout = [varargout {axon_data}];
  else varargout = [varargout {axon_file}]; 
  end
end

if any(request('fascicle')) % just the fascicle requested, not finalised 
    
  if exist('axon_data','var'), F_is_embedded = isfield(axon_data,'nerve'); 
  elseif exist('axon_file','var')  
    file_info = whos('-file',axon_file); 
    F_is_embedded = any(strcmp({file_info.name},'nerve'));
  elseif any(named('-anat')) % user supplied    
    F_is_embedded = true; 
    varargout = [varargout {get_('-anat')}];
  else    
    error('TODO: requested anatomy but not sure how to satisfy request')
  end

  if ~F_is_embedded % older axons files might not have this embedded
    f_list = dir(tools.file('~/source/fascicles/*.splines.dat'));
    nerve = mesh.read_dat_file(tools.INPUT_file(f_list, p_(axon_file)));
    
    if exist('axon_data','var'), varargout{end}.nerve = nerve;
    else
      error('TODO: not sure what the best option is here. %s', ... 
            'Probably easier to just run models.axon_population again')
    end
  end
end

%% Find /stimulation/ folder (if requested)
if any(request('stim'))

  if any(named('-stim'))
    stim_folder = get_('-stim');
    if any(stim_folder == '~'), stim_folder = tools.file(stim_folder); end
    if exist(stim_folder,'dir') ~= 7
      error('stim folder %s not found', stim_folder)
    end
  else
    s_list = dir(tools.file('stim~/'));
    s_list(cellfun(@(n) all(n=='.'),{s_list.name})) = []; 
    
    if exist('eidors_file','var')
      [~,e_file] = fileparts(eidors_file);
    else e_file = par_('stim');      
      disp('TODO: load /stimulation without /eidors')
    end
    
    n_part = regexp(e_file,'\([^\)]*\)','match','once'); 

    if isempty(s_list)
      error('Please run models.nerve_stimulation')
    end

    if any(strcmp({s_list.name}, e_file))
      sel = strcmp({s_list.name}, e_file);
    elseif any(contains({s_list.name}, n_part))
      sel = contains({s_list.name}, n_part);
    else
      sel = ([s_list.datenum] == max([s_list.datenum])); 
      warning('ViNERS:stimFileName', '%s not found in %s, using %s', ...
              n_part, s_list(1).folder, s_list(sel).name)
    end  
    if sum(sel) > 1 % get newest
      sel = ([s_list.datenum] == max([s_list(sel).datenum])); 
    end

    stim_folder = p_(s_list(sel)); 
  end
  
  varargout = [varargout {stim_folder}];
end


%% Find /wave/ folder (if requested)
if any(request('wave'))

  wave_folder = tools.file('waves~\');
  
  % plot_functions accept wave_folder as a first argument
  if numel(user_args) > 0, wave_folder = user_args{1}; end
  % plot_functions accept just the name of the folder inside of \waves\
  if ~any(ismember(wave_folder,'/\')), ...
      wave_folder = ['waves~\' wave_folder]; 
  end
  
  % models.function( ..., '-waves', eidors_file ) also valid unless excluded
  if any(named('-waves')) && ~any(request('w-positional-only'))
         wave_folder = get_('-waves'); 
  end
  if any(wave_folder == '~'), wave_folder = tools.file(wave_folder); end 
  
  if strncmp(wave_folder,'?',1) % user requested a menu
    
    if numel(wave_folder) == 1, wave_folder = ['?' par_('wave')]; end % '?'
    wave_folder = uigetdir(wave_folder(2:end));
    
  elseif exist(wave_folder,'dir') ~= 7 || ... % use handle or menu anyway
         strcmp(wave_folder,tools.file('waves~\')) % just the outer folder?
    
    wave_folder = par_('wave');    
    if isa(wave_folder,'function_handle')
      % e.g. 'eidors',@() tools.file('get','eidors~/stim*.mat','newest')
      wave_folder = wave_folder(); % user_args{:}
    else
      if any(wave_folder == '~'), wave_folder = tools.file(wave_folder); end
      wave_folder = uigetdir(wave_folder);
    end
  end
  
  if opts.do_multi, error wave_multiselect, end
  varargout = [varargout {wave_folder}];
  
  if any(request('wave+context')) % uses tools.INPUT_file
    
    if ~any(request('eidors')) % already requested
      
      f = dir(tools.file('sub~\eidors\sensitivity*.mat'));
      if any(named('-ei')), eidors_file = get_('-ei');
      else eidors_file = tools.INPUT_file(f, wave_folder);  
      end
      
      if opts.do_LOAD, announce_load_(eidors_file)
        EM = load(eidors_file); [~,EM.filename] = fileparts(eidors_file);
        varargout = [varargout {EM}];   
      else varargout = [varargout {eidors_file}];   
      end
    end
    if ~any(request('axon')) % already requested
      
      axons_folder = tools.file('axons~\');
            
      if any(named('-ax')), 
        axon_file = get_('-ax');
        if any(axon_file == '?'), 
          [af.name,af.folder] = uigetfile('axons*.mat',[],tools.file('axons~'));     
           axon_file = p_(af);
        end

      else % use defaults (tools.INPUT_file)
        axon_file = dir([axons_folder filesep '*.mat']);  
        if isempty(axon_file), error('please run models.axon_poplation'), end
        axon_file = tools.INPUT_file(axon_file, eidors_file);
      end

      if opts.do_LOAD, announce_load_(axon_file)
         EM.axons = load(axon_file);
         varargout{end} = EM; % overwrite, assuming position
         
         if isfield(EM.axons,'nerve')
           nerve = EM.axons.nerve;
         else
           f = dir(tools.file('~/source/fascicles/*.splines.dat')); 
           spline_file = tools.INPUT_file(f, axons_file);
           nerve = read_dat_file(spline_file,'index',EM.info.SplineIndex);           
         end
         
         if any(named('-delta-xy'))
           error TODO_delta_XY
%             dxy = get_('-delta-xy');
%             fprintf('Displacing fascicle by [%0.3f, %0.3f] mm\n', dxy)
%             axon_data.pop.fibre_xy = axon_data.pop.fibre_xy + dxy;
%             axon_data.pop.unmyelinated_xy = axon_data.pop.unmyelinated_xy + dxy;
%             
%             nerve.fascicle = nerve.outline / 1000; % convert to mm
%             nerve.fascicle(:,2,:) = nerve.fascicle(:,2,:) + 0.050; % add offset
         end
           varargout = [varargout {nerve}];
      else varargout = [varargout {axon_file}];
           f = dir(tools.file('~/source/fascicles/*.splines.dat')); 
           spline_file = tools.INPUT_file(f, axon_file);
           varargout = [varargout {spline_file}];
      end
    end
  end  
end





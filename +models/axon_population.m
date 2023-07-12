
function [pop, sam] = axon_population(varargin)
% pop = make_axon_population([dataset_name], ... )
% make_axon_population takes a set of fascicle contours and a specification
%   for the axon populations, and populates the fascicles with axons with
%   varying diameters and g-ratios. 
% 
% We have data-sets scraped from published papers or the SPARC data portal. 
%  Data-sets scraped from published papers often rely on reconstructing the 
%  underlying joint distribution from the published marginal histograms.
%  For a list of valid data-sets, run models.axon_population('help data')
% 
% by default, the current fascicle pattern defined in 
% >> mesh.insert_gmsh_fascicles
%   is used to define the fascicle, and `axonType.csv` is read to define 
%   the total axon count (default) or axon count per mm2 (if axonType.csv
%   is defined in terms of "per_mm2"). 
% 
% EXAMPLE:
%
% mesh.insert_gmsh_fascicles('-setup','-file','nerve.xml', ...
%                            '-config',options)
% models.make_axon_population('pelvic nerve')
% 
% 
% BASIC OPTIONS:
%  -output 'file' : save the resulting axon population to the specified
%                   output file.
%  -branch 'name' : some axon populations have more detail provided for 
%                   individual nerve branches. This selects which branch to
%                   utilise, otherwise the average composition of the 
%                   branches will be used.
%  -case 'name'   : as -branch
%  -upsample [amount] : up-sample axon counts by [amount] e.g. '-ups',1.2 
%                       gives 20% more axons in generated file
%  -down [amount] : down-sample axons by [amount]. Values greater than 1 
%                   are interpreted as factors e.g. -down 2 and -down 0.5 
%                   result in the same outcome: half as many axons. 
%  -file [file]   : define source axons file for axon subtype definitions.
%                   Otherwise, new clusters are generated and you will be 
%                   required to run models.membrane_currents
%  -ngroups [24]  : define # of axon groups for sub-type clustering to
%                   accellerate computation of membrane currents and SFAPs.
%  -class [table] : define # of axons of each class. By default, the 
%                   following table is used: 
% 
%  t.Type = {'Myelinated Afferent'; 'Myelinated Efferent';
%            'Unmyelinated Afferent'; 'Unmyelinated Efferent'}
%  t.Count = [475;       395;    1200;   2800]
%  t.Model = {'Gaines'; 'MRG'; 'Sundt'; 'Sundt'}
%  (based on Hulsebosch doi://10.1002/cne.902110102)
% 
% For a complete list of options (including options provided for historical
%   compatability), run models.axon_population('help options').
% 
% For a complete list of the available data, run 
%   models.axon_population('help data')
% 
% For more information about adding new datasets to ViNERS, see <link>
%
% Version 0.3, refactored to enable code-free addition of new datasets 
% Version 0.2, refactored to merge approaches
% Calvin Eiber 16-April-2020 ceiber@unimelb.edu.au


close all
if nargin == 0, varargin = {'SN','!!!do-types-warning'}; end

varargin = tools.opts_to_args(varargin,'axons');
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('help')) || any(named('-h')), display_help_menu(named), return, end

%% Get source data specification and nerve outline to fill in

[D,pop,nom,axon_types] = gather_sourcedata(varargin{:});

nerve = get_nerve_outline(varargin{:}); 
is_column = @(t,v) any(strcmp(t.Properties.VariableNames,v));

if any(named('-tiger'))
    axon_types.Model(strcmp(axon_types.Model,'Sundt')) = {'Tigerholm'};
end

%% How many axons? (note: number of axons can be downsampled later)
if is_column(axon_types,'per_mm2')   
  if ~(any(named('-use-count')) && is_column(axon_types,'Count'))
    f_area_mm2 = arrayfun(@(f) polyarea(nerve.fascicles(:,1,f), ...
                                        nerve.fascicles(:,2,f)), ...
                                 1:size(nerve.fascicles,3));
    axon_types.Count = round(axon_types.per_mm2  * sum(f_area_mm2));
  end
end

if any(named('-ups')), scale_factor = get_('-ups'); % Up-sample axon_types
  % scale_factor = 10.^(abs(log10(scale_factor)));
  axon_types.Count = round(axon_types.Count .* reshape(scale_factor,[],1)); 
end
if any(named('-down')), scale_factor = get_('-down'); % Down-sample axon_types  
  scale_factor = 10.^(-abs(log10(scale_factor)));
  axon_types.Count = round(axon_types.Count .* reshape(scale_factor,[],1)); 
  axon_types.Count(axon_types.Count < 1) = 1;   
  warning('ViNERS:make_axon_population:earlyDownsample',...
           'Downsampling axon population by %g, %s. %s %s', 1./scale_factor(1), ...
           'computed recordings or ECAPs may not be accurate', ...
           'Downsampling can be accomplished in later modules which', ...
           'better preserves whole-nerve estimates of signal amplitudes.')
end
clear scale_factor

%%
if ~isempty(pop) % Construct population from histogram data   
  % We have the underlying EM data-set, resample so it has the correct # 
  % of axons given our chosen fascicles
  pop = sample_from_population(D,pop,axon_types);
else 
  % Make populations from histograms 
  pop = generate_from_histograms(D,nom,axon_types); 
end

if ~any(named('-no-fig')), figure(1), mk_hist_scatter_figure(D,pop); end

if nargout == 1, return, end % Load or generate /just/ the population part of the data 

%% Down-sample the data 

nG = 24; % <<< Control down-sample resolution 
if any(named('-nG')), nG = varargin{find(named('-nG'))+1}; end

[pop,sam] = arrange_data_sample(pop, nG);

color = summer(nG); % This is the colorscheme we're adopting for these
color(2:2:end,:) = 1-color(2:2:end,:); 
color = color([2 11 10 1],:).*[.8 .8 1 1]';

nF = size(nerve.fascicles,3);

if ~any(named('-no-fig'))
  %%
  figure(2), clf, hold on, % C = lines(nF);
  set(gcf,'Position',get(gcf,'Position') .* [1 .7 1.3 1])
  if iscell(nerve.outline)
    fill(nerve.outline{end}(:,1),nerve.outline{end}(:,2),[.9 .9 .9],'EdgeColor','none')
  end  
  for ff = 1:nF
      fill(nerve.fascicles(:,1,ff),nerve.fascicles(:,2,ff),'w','EdgeColor',[.3 .3 .3],'LineWidth',1.2)    

      sel = (f_id == ff) & ~is_myelin & is_afferent;
      plot(xy(sel,1),xy(sel,2),'.','Color',color(1,:),'MarkerSize',10) 

      sel = (f_id == ff) & ~is_myelin & ~is_afferent;
      plot(xy(sel,1),xy(sel,2),'.','Color',color(2,:),'MarkerSize',10) 

      sel = (f_id == ff) & is_myelin & is_afferent;
      plot(xy(sel,1),xy(sel,2),'o','MarkerFaceColor',(color(3,:)+1.5)/2.5, ...
                  'MarkerSize',5,'Color',color(3,:),'LineWidth',1.2)

      sel = (f_id == ff) & is_myelin & ~is_afferent;
      plot(xy(sel,1),xy(sel,2),'o','MarkerFaceColor',(color(4,:)+1.5)/2.5, ...
                  'MarkerSize',5,'Color',color(4,:),'LineWidth',1.2)
  end
  axis image, tools.tidyPlot

  l_str = {'Fascicle outline','unmyelinated afferent', ...
             'unmylinated efferent','myelinated afferent', ...
                                    'myelinated efferent'};
  for tt = height(axon_types):-1:1    
    if axon_types.Count(tt) == 0, l_str(tt+1) = []; end
  end
  
  if iscell(nerve.outline), l_str = [{'Endoneurium'} l_str]; end
  legend(l_str{:},'location','best')
end

%% clean up and save

assert(all(sam.A_diam > 0),'in downsample, some A axons had diameter <= 0. Try reducing -nG')
assert(all(sam.C_diam > 0),'in downsample, some C axons had diameter <= 0. Try reducing -nG')
assert(all(sam.A_inner > 0),'in downsample, some A axons had g-ratio <= 0. Try reducing -nG')

clear ok in_a_bv vp perimeter_fcn ad_r bb c2 cat cf cp est_fit
clear est_init fd_r gr gr_ok gr_r graph_roi lhs lhs_AD lhs_FD lhs_GR n nr
clear lloyd_relax N n_ratio pp rhs rhs_out x y aidx v1 v2 v3 v4 v5

if max(abs(pop.fibre_xy(:))) > 10 && ~any(named('-units-no-c'))
  pop.fibre_xy = pop.fibre_xy / 1e3;
  pop.unmyelinated_xy = pop.unmyelinated_xy / 1e3;
    warning('ViNERS:make_axon_population:outputUnits', ...
            'Converting output units from µm to mm.')
end

if any(named('-mat-legacy'))
  
  if mod(height(axon_types),2)
       axon_types = reshape([axon_types.Count;0],2,[]); 
  else axon_types = reshape( axon_types.Count, 2, []); 
  end  
     varnamelist = {'pop','sam','nerve','opts'};
else varnamelist = {'pop','nerve','opts'};
end

pop.axon_table = axon_types;
pop.axon_color = color;

pop.units.axon_xy = 'mm'; 
pop.units.axon_diameter = 'um (inner diameter)'; 
pop.units.fibre_diameter = 'um (outer diameter)'; 
pop.units.g_ratio = '(ratio of diameters)'; 
pop.units.fascicle_index = ['nF = 1..' num2str(nF)];
pop.units.axon_table = 'total axon count';

opts = varargin;  %#ok<NASGU> 

if any(named('-out')), output = get_('-out'); 
% elseif any(named('-source')), output = get_('-file'); 
else output = tools.file('get','axons~/axons (%s).mat','next',D.name);
end
if ~any(ismember(output,'/\()'))
    output = ['sub~/axons/axons (' output ').mat'];
end

if any(output == '~'), output = tools.file(output); end
if ~any(named('-mat-legacy')) || any(named('-csv')) || any(named('-xls'))
  pop = convert_output_structure(pop,sam);
end

if any(named('-no-save')), return, end

if any(named('-csv')), export_tables(pop,sam,nerve,'csv',output)  
elseif any(named('-xls')), export_tables(pop,sam,nerve,'xlsx',output)
else  
  if ~isfolder(fileparts(output)), mkdir(fileparts(output)); end
  while exist(output,'file'), output = strrep(output,'.mat','_NEW.mat'); end
  fprintf('Saving %s\n',tools.file('T',output))
  save(output, varnamelist{:}); 
end

%%
return

%% HELP menu


function display_help_menu(named)

d = @disp;   
if any(named('options')) || any(named('help o'))

d(' ')  
d('pop = make_axon_population([type], ... )')
d(' ')
d('BASIC OPTIONS:')
d('-output ''file'' : save the resulting axon population to the specified')
d('                 output file.')
d('-anat [nerve]  : Use the specified nerve anatomy (where "nerve" is')
d('                 the output of plots.preview_fascicles)')
d('-branch ''name'' : some axon populations have more detail provided for')
d('                 individual nerve branches. This selects which branch to')
d('                 utilise, otherwise the average composition of the') 
d('                 branches will be used.')
d('-case ''name''   : as -branch')
d('-upsample [amount] : up-sample axon counts by [amount] e.g. ''-ups'',1.2')
d('                     gives 20% more axons in generated file')
d('-down [amount] : down-sample axons by [amount]. Values greater than 1 ')
d('                 are interpreted as factors e.g. -down 2 and -down 0.5 ')
d('                 result in the same outcome: half as many axons. ')
d('-remake [file] : define source axons file for axon subtype definitions.')
d('                 Otherwise, new clusters are generated and you will be ')
d('                 required to run models.membrane_currents')
d('-ngroups [24]  : define # of axon groups for sub-type clustering to')
d('                 accellerate computation of membrane currents and SFAPs.')
d('-class [table] : define # of axons of each class. By default, a table')
d('                 based on Hulsebosch doi://10.1002/cne.902110102 is used.')
d(' ')
d('SAVE OPTIONS:')
d('-csv        : Generate a set of .CSV tables instead of the usual .mat')
d('              file output. These will not be usable with ViNERS but may')
d('              be useful for other scientific purposes.')
d('-xlsx       : as .CSV, but generate .XLSX tables instead.')
d('-mat-legacy : Generate a .mat file with old variable names and data')
d('              structures (provided for backwards-compatibility with')
d('              earlier versions of ViNERS).')
d('-no-save    : Do not save the output file (preview mode).')
d(' ')
d('ADVANCED OPTIONS:')
d('-axon-path [path] : If set, read axon populations from a location')
d('                    other than the default (~\source\axons\)')
d('-backup [file] : If not enough data is available for the selected nerve,')
d('                 read additional histogram data from the specified file.')
d('-indexed-only  : If set, only axon populations which are indexed in')
d('                 ~\source\axons\index.json can be used. ')
d('                 (default: all folders in ~\source\axons\ can be used)')
d('-keep          : If set, override number of axons generated')
d('-n-Lloyd [n]   : Set the number of Lloyd relaxation steps. We actually')
d('                 found that the axon layout resembled our nerve sections')
d('                 best with n = 0, which is the default value')
d('-q             : If set, quiet mode')
d('-secondary []  : synonym for -backup')
d('-units-no-convert   : If not set, if any axons have |x| or |y| > 10 it')
d('                      is assumed that the input geometry was specified')
d('                      in µm and length units are converted to mm.')
d('-use-count : If set, absolute counts for each axon type will be used to')
d('             determine the number of axons if available (default: use')
d('             per_mm2 values if available).')
d('-zscore    : cluster axons by (z-scored) fibre diameter and g-ratio.')
d('             Default: cluster by fibre and axon diameter. We found that')
d('             clustering on fibre and axon diameter converged sooner')
d('             (in the sense of changes in simulated ECAPs) and so this')
d('             is not recommended, but provided for backwards compatibility.')
d(' ')
d('Version 1.1, 2021 November 16. Upgrades to enable membrane current reuse')
d('             Calvin Eiber      ceiber@unimelb.edu.au')
d(' ')

elseif any(named('data')) || any(named('help d'))
  %% 
  
  index = tools.parse_json(tools.file('~/source/axons/index.json'));
  index = [index.list{:}];

  folder_names = dir(tools.file('~/source/axons/*-*'));
  
  if ~any(named('-indexed-only'))
    % append literal names for folders not in index.json    
    folder_names(ismember({folder_names.name},{index.entry})) = []; 
    append = index(ones(size(folder_names))); 
    [append.short] = deal(folder_names.name);
    [append.entry] = deal(folder_names.name);
    index = [index append]; 
  end

  [uie,~,synonym] = unique({index.entry}); 
  
  
  d(' ')
  d('Axon distributions available for simulation:')
  d(' ')
  fprintf(' %-20s\t: Count : My Axons   : Unmy axons\n', 'Nerve')
  d('----------------------------------------------------------------')
  
  p_ = @(x,varargin) [x.folder filesep x.name varargin{:}];  
  
  for ii = 1:numel(uie)
    
    list = dir([folder_names(1).folder filesep uie{ii}]);
    
    has_2DM = any(contains({list.name},'MyelinatedAxons.csv'));
    has_UD  = any(contains({list.name},'UnmyelinatedAxons.csv'));
    has_AT  = any(strcmp({list.name},'axonType.csv')); 
    has_H   =    (strncmp({list.name},'histograms',8)); 
    
    has_GR = has_2DM; 
    has_MW = has_2DM; 
    
    if ~has_2DM && any(has_H)
      
      t = readtable(p_(list(find(has_H,1)))); 
      n = lower(t.Properties.VariableNames); 
      
      if ~has_UD, has_UD = any(strncmp(n,'unmyelinated',8)); end      
      if ~has_GR, has_GR = any(contains(n,'ratio')); end
      if ~has_MW, has_MW = any(contains(n,'width')); end
      
      clear t n
    end
    
    if has_AT, sc = 'true '; else sc = 'false'; end
    
    if has_2DM,    sm = '[3] D vs G';
    elseif has_GR, sm = '[2] D, G'; 
    elseif has_MW, sm = '[2] D, W'; 
    else           sm = '[1] D only'; 
    end
    
    if has_UD, su = 'true '; else su = 'false'; end
    
    fprintf(' %-20s\t: %5s : %-10s : %5s \n', uie{ii}, sc,sm,su)
    
    if sum(synonym == ii) > 1
      sy = unique({index(synonym == ii).short});
      sy(strcmp(sy,uie{ii})) = [];
      sy = sprintf(', %s', sy{:});
      fprintf('      [alias: %s]\n', sy(2:end));
    end
    
  end % for each folder 
else
  
  help('models.axon_population')
end
  
  
  
  

return



%% parse input arguments to get source data
function [D,pop,nom,tab] = gather_sourcedata(varargin)

pop = []; 

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

index = tools.parse_json(tools.file('~/source/axons/index.json'));
index = [index.list{:}];

if ~any(named('-indexed-only'))
  % append literal names for folders not in index.json
  folder_names = dir(tools.file('~/source/axons/*-*'));
  folder_names(ismember({folder_names.name},{index.entry})) = []; 
  append = index(ones(size(folder_names))); 
  [append.short] = deal(folder_names.name);
  [append.entry] = deal(folder_names.name);
  index = [index append]; 
end

if ~any(named('-q')), print = @fprintf; else print = @(varargin) []; end

if any(named('-axon-path')), source_path = get_('-axon-path');
else 
  %% Parse input arguments 

  arg = upper(varargin(cellfun(@ischar,varargin)));

  sel = 0; 

  for vec = [ {index.entry}, {index.short}, strcat('-',{index.short})] % fixed CE 13-07-23 for R2022a

    hit = ismember(arg,upper(vec{1}));
    if any(hit), sel = find(hit,1); break, end
  end

  if ~any(sel) || any(named('!!!do-types-warning'))
    [e,x] = unique({index.entry});
    for ii = 1:numel(x)
      fprintf(2,'-%-5s : %s\n', index(x(ii)).short,e{ii})
    end    
    if ~any(sel) 
      error('none of the above sets of installed axon parameters selected')
    end
  end

  sel = find(ismember(upper(vec{1}),arg{sel}),1);
  
  source_path = tools.file(['~/source/axons/' index(sel).entry]);  
end

print('Loading axons from %s\n', tools.file('T',source_path));


list = dir(source_path);
list([list.isdir]) = []; 

D = struct;
D.name = index(sel).entry;
nom = {}; 

%% Import axon tables if available 

p_ = @(x) [x.folder filesep x.name];
MY_csv_files = contains({list.name},'MyelinatedAxons.'); 

if any(MY_csv_files)
  
  t = arrayfun(@(x) readtable(p_(x)), list(MY_csv_files),'unif',0); 
  pop = struct; 
  
  fibre_props = {'fibre_diam','fibre_axon','fibre_gratio','fibre_xy'};
  for v = fibre_props,pop.(v{1}) = []; end
  
  for ff = 1:numel(t) % import each table 
    
    done = false(1,3); % diam axon gratio xy

    var = lower(t{ff}.Properties.VariableNames);
    if any(contains(var,'minferet'))
      minFeret = contains(var,'minferet');
      
      sel = contains(var,'fib') & contains(var,'diam') & minFeret;
      if any(sel)
        if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for fibre_diameter_minFeret, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(1) = true;
        pop.fibre_diam = [pop.fibre_diam; t{ff}{:,sel}];
      end

      sel = contains(var,'ax') & contains(var,'diam') & minFeret;
      if any(sel)      
        if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for axon_diameter_minFeret, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(2) = true;
        pop.fibre_axon = [pop.fibre_axon; t{ff}{:,sel}];
      end

      sel = contains(var,'ratio') & minFeret;
      if any(sel)
        if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for g_ratio_minFeret, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(3) = true;
        pop.fibre_gratio = [pop.fibre_gratio; t{ff}{:,sel}];
      end
    end % minFeret
    
    if ~all(done)
      sel = contains(var,'fib') & contains(var,'diam');
      if any(sel) && ~done(1)
        if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for fibre_diameter, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(1) = true;
        pop.fibre_diam = [pop.fibre_diam; t{ff}{:,sel}];
      end

      sel = contains(var,'ax') & contains(var,'diam');
      if any(sel) && ~done(2)      
        if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for axon_diameter, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(2) = true;
        pop.fibre_axon = [pop.fibre_axon; t{ff}{:,sel}];
      end

      sel = contains(var,'ratio');
      if any(sel) && ~done(3)
        if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for g_ratio, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(3) = true;
        pop.fibre_gratio = [pop.fibre_gratio; t{ff}{:,sel}];
      end
    end % generic    
    
    sel = contains(var,'_x') | strcmp(var,'x'); 
    if any(sel)
      if sum(sel) > 1
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for X, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
      end
      sel = find(sel,1) + [0 1]; % assert
      
      if ~strcmp(var{sel(2)},strrep(var{sel(1)},'x','y'))
          warning('ViNERS:make_axon_population:XYorder', ...
                  'the column after %s was %s (expected %s)', ...
                   t{ff}.Properties.VariableNames{sel}, strrep(var{sel(1)},'x','y'))
      end
      
      pop.fibre_xy = [pop.fibre_xy; t{ff}{:,sel}];
    end
    
    if ~all(done)
      
      my_list = list(MY_csv_files);
      
      warning('ViNERS:make_axon_population:missingData', ...
              'a required data column was missing from %s, skipping ...', ...
                   my_list(ff).name)
      n = min(cellfun(@(v) size(pop.(v),1),fibre_props));
      for v = fibre_props'
        pop.(v{1}) = pop.(v{1})(1:n,:);
      end
    end
  end
  
  if size(pop.fibre_xy,1) ~= size(pop.fibre_axon,1)
    pop.fibre_xy = []; % throw out without comment
  end
  
  if isempty(pop.fibre_axon), pop = []; end
  nom = [nom; {'pop.fibre_diam','fibre_diam'; ...
               'pop.fibre_axon','fibre_axon'; ...
               'pop.fibre_gratio','fibre_gratio'}];
end

UNMY_csv_files = contains({list.name},'UnmyelinatedAxons.'); 

if any(UNMY_csv_files)
  
  t = arrayfun(@(x) readtable(p_(x)), list(UNMY_csv_files),'unif',0); 
  if isempty(pop), pop = struct; end
  
  fibre_props = {'unmyelinated_diam','unmyelinated_xy'};
  for v = fibre_props,pop.(v{1}) = []; end
  
  for ff = 1:numel(t) % import each table 
    
    done = false(1,1); % diam axon gratio xy
    var = lower(t{ff}.Properties.VariableNames);
    
    if any(contains(var,'minferet'))
      minFeret = contains(var,'minferet');
      
      sel = contains(var,'diam') & minFeret;
      if any(sel),
        if sum(sel) > 1,
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for fibre_diameter_minFeret, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(1) = true;
        pop.unmyelinated_diam = [pop.unmyelinated_diam; t{ff}{:,sel}];
      end
    end % minFeret
    
    if ~all(done)
      sel = contains(var,'diam');
      if any(sel) && ~done(1)
        if sum(sel) > 1,
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for fibre_diameter, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
        end
        sel = find(sel,1); done(1) = true;
        pop.fibre_diam = [pop.fibre_diam; t{ff}{:,sel}];
      end
    end % generic    
    
    sel = contains(var,'x') & ~contains(var,'ax');
    if any(sel)
      if sum(sel) > 1, 
          warning('ViNERS:make_axon_population:multipleDataColumns', ...
                  '%d entries found for X, using %s', ...
                   sum(sel),t{ff}.Properties.VariableNames{find(sel,1)})
      end
      sel = find(sel,1) + [0 1]; % assert
      
      if ~strcmp(var{sel(2)},strrep(var{sel(1)},'x','y')), 
          warning('ViNERS:make_axon_population:XYorder', ...
                  'the column after %s was %s (expected %s)', ...
                   t{ff}.Properties.VariableNames{sel}, strrep(var{sel(1)},'x','y'))
      end
      
      pop.unmyelinated_xy = [pop.unmyelinated_xy; t{ff}{:,sel}];
    end
    
    if ~all(done)
      
      my_list = list(UNMY_csv_files);
      
      warning('ViNERS:make_axon_population:missingData', ...
              'a required data column was missing from %s, skipping ...', ...
                   my_list(ff).name)
      n = min(cellfun(@(v) size(pop.(v),1),fibre_props));
      for v = fibre_props'
        pop.(v{1}) = pop.(v{1})(1:n,:);
      end
    end
  end
  
  if size(pop.unmyelinated_xy,1) ~= size(pop.unmyelinated_diam,1)
    pop.unmyelinated_xy = []; % throw out without comment
  end

  if isempty(pop.unmyelinated_diam)
    pop.unmyelinated_diam = []; 
    pop.unmyelinated_xy = [];  
  else nom = [nom; {'pop.unmyelinated_diam','unmyelinated_diam'}];  
  end
elseif any(MY_csv_files)
  pop.unmyelinated_diam = []; 
  pop.unmyelinated_xy = [];
end

if any(named('-case')),       branch_ID = get_('-case');
elseif any(named('-branch')), branch_ID = get_('-branch');
else                          branch_ID = ''; 
end
  

%% Import histograms flexibly from histograms.csv
if isempty(pop) || isempty(pop.fibre_axon) || isempty(pop.unmyelinated_diam) % else 
  sel = contains({list.name},'histograms');  
  
  if ~any(sel)
    if isempty(pop) || isempty(pop.fibre_axon)
      
        missingCSV = 'MyelinatedAxons.csv';
        error('if %s are not supplied, histograms.csv must be supplied %s', ... 
                missingCSV, '(note: this is case-sensitive).')

    end
    
    if any(named('-sec')), supp_file = get_('-sec'); 
    elseif any(named('-backu')), supp_file = get_('-backu'); 
    else
      supp_file = tools.file('~/source/axons/rat-vagus-cervical/histograms.csv');
      warning('ViNERS:make_axon_population:missingHistogramCSV',...
              'histograms.csv missing, supplementing data with %s', ...
              supp_file)
    end
    sel(end+1) = true;
   [list(sel).folder,list(sel).name, sfe] = fileparts(supp_file);
    list(sel).name = [list(sel).name sfe];
  end
  
  if sum(sel) > 1
   
    sel = sel & contains(lower({list.name}),lower(branch_ID));
    if any(sel), 
      sel(sel) = ([list(sel).bytes] == max([list(sel).bytes]));
    end
    
    if any(named('-case')) && any(named('-branch'))
         branch_ID = get_('-branch');
    else branch_ID = ''; 
    end
    
  end
  
  print('Loading %s\n',tools.file('T',p_(list(find(sel,1)))))
  hfile = readtable(p_(list(find(sel,1))));
    
  done = false(size(hfile.Properties.VariableNames));
  
  for cc = 1:numel(done)
    %%
    var = lower(hfile.Properties.VariableNames{cc});
    if done(cc), continue, end
    if ~isnumeric(hfile{:,cc}), continue, end    
    ok = ~isnan(hfile{:,cc});
    
    D.(var).x = hfile{:,cc}(ok);
    D.(var).y = hfile{:,cc+1}(ok);
    
    % normalise
    D.(var).y = D.(var).y-min([0;D.(var).y]);
    D.(var).y = D.(var).y / sum(D.(var).y);
    
    done(cc+[0 1]) = true; 
    
    if contains(var,'axondiam') && contains(var,'unmy')
                                      nom = [nom;{ var 'unmyelinated_diam'}]; %#ok<AGROW>
    elseif contains(var,'fibrediam'), nom = [nom;{ var 'fibre_diam'}];        %#ok<AGROW>
    elseif contains(var,'fiberdiam'), nom = [nom;{ var 'fibre_diam'}];        %#ok<AGROW>
    elseif contains(var,'axondiam'),  nom = [nom;{ var 'fibre_axon'}];        %#ok<AGROW>
    elseif contains(var,'gratio'),    nom = [nom;{ var 'fibre_gratio'}];      %#ok<AGROW>
    elseif contains(var,'g_ratio'),   nom = [nom;{ var 'fibre_gratio'}];      %#ok<AGROW>
    elseif contains(var,'width'),     nom = [nom;{ var 'fibre_mywidth'}];      %#ok<AGROW>      
    elseif cc < width(hfile)
      if ~contains(lower(hfile.Properties.VariableNames{cc+1}),'count')
        done(cc+[0 1]) = false;
      end
      warning('ViNERS:make_axon_population:histogramUnknown',...
              'Not sure what to do with column "%s", skipping ...', ...
              hfile.Properties.VariableNames{cc})
      continue
    end
  end % loop over columns
  
  if ~isempty(branch_ID) % Get (internal) sub-sample for branch or case ID 
        
    sel = contains(lower(nom(:,1)),lower(['_' branch_ID]));     
    if ~any(sel)
      warning('ViNERS:make_axon_population:branchUnknown', ...
              'Nothing in %s matched "%s", ignoring...', ...
               tools.file('T',source_path), branch_ID)
    else
      sel = sel | ismember(nom(:,2), setdiff(nom(:,2),nom(sel,2)));
      sel = sel | contains(nom(:,1),'pop.');
      nom = nom(sel,:); 
    end
  end
  
  if any(strncmp(nom(:),'aff',3))
    [D,nom] = parse_affeff_balance(D,nom); 
  end
   
  %% Ensure enough data to produce a population has been loaded 
  
  n_supplied = @(n) sum(strncmp(unique(nom(:,2)),n,numel(n)));  
  is_adequate = n_supplied('fib') > 1 & n_supplied('unmy') > 0;
  if ~is_adequate
    
    missing_field = {}; 
    if n_supplied('unmy') <= 0
      missing_field{end+1} = 'unmyelinated_diam';
    end    
    if n_supplied('fib') <= 1
      missing_field{end+1} = 'fibre_gratio';
    end
    assert(~isempty(missing_field),'unclear what is missing')    
    
    if any(named('-sec')), supp_file = get_('-sec'); 
    elseif any(named('-backu')), supp_file = get_('-backu'); 
    else
      supp_file = tools.file('~/source/axons/rat-vagus-cervical/histograms.csv');
      warning('ViNERS:make_axon_population:supplementingData',...
              '%s not set, loading a histogram at random from %s', missing_field{1}, supp_file)
    end
    
    supp_file = readtable(supp_file);
    
    if any(strcmp(missing_field,'fibre_gratio'))
      ok = ~isnan(supp_file.Myelinated_gRatio);
      D.g_ratio.x = supp_file.Myelinated_gRatio(ok);
      D.g_ratio.y = supp_file.Myelinated_gRatioCount(ok);
      nom = [nom; {'g_ratio', 'fibre_gratio'}];
    end
    if any(strcmp(missing_field,'unmyelinated_diam'))      
      ok = ~isnan(supp_file.Unmyelinated_axonDiam);
      D.unmyelinated_axondiam.x = supp_file.Unmyelinated_axonDiam(ok);
      D.unmyelinated_axondiam.y = supp_file.Unmyelinated_axonCount(ok);
      nom = [nom; {'unmyelinated_axondiam', 'unmyelinated_diam'}];
    end
  end
 
end % histogram load code

%% Load axon populations table
sel = contains({list.name},'axonType.');  
if any(named('-class-pn')) || any(named('-class-default'))
  % Afferent; Efferent (A-delta C-fibre), from H+C Table 6
  % DOI: 10.1002/cne.902110102
  tab = [475  1200;  275+120   2100+700]; 
elseif any(named('-class')), tab = get_('-class');
elseif ~any(sel) 
  if ~isempty(pop)
       error('TODO: implement types table from imported population')
  else error('if MyelinatedAxons.csv are not supplied, histograms.csv must be supplied.')  
  end
else tab = readtable(p_(list(find(sel,1))));  
end

if isnumeric(tab)  
  if ~all(size(tab) == 2)
    warning('ViNERS:make_axon_population:tableSize', ...
            'Expected 2x2 table: [A-type aff, C-type aff; A-type eff, C-type eff]')
  end
  
  tab = round(tab); 
  
  t = struct;
  t.Type = repmat({'unknown type'},[numel(tab) 1]);
  t.Count = tab(:);
  t.Model = repmat({''},[numel(tab) 1]);
  if numel(tab) >= 1, t.Type{1} = 'Myelinated Afferent'; 
                      t.Model{1} = 'Gaines';              end
  if numel(tab) >= 2, t.Type{2} = 'Myelinated Efferent'; 
                      t.Model{2} = 'MRG';                 end
  if numel(tab) >= 3, t.Type{3} = 'Unmyelinated Afferent'; 
                      t.Model{3} = 'Sundt';               end
  if numel(tab) >= 4, t.Type{4} = 'Unmyelinated Efferent'; 
                      t.Model{4} = 'Sundt';               end
  tab = struct2table(t);
end


v  = fieldnames(D);
if isempty(v)
  D.g_ratio.x = (1:10)/10; % noted by mk_hist_scatter_figure
  v  = fieldnames(D);
end

% Parse the "D" object to determine how we're going to generate the sample
sel = contains(v,'ratio');
if any(sel), D.index(1) = find(sel,1); else D.index(1) = -1; end % g ratio
sel = contains(v,'diam') & contains(v,'fi') & ~contains(v,'unmy');
if any(sel), D.index(2) = find(sel,1); else D.index(2) = -1; end % myel fibre diam
sel = contains(v,'diam') & contains(v,'ax') & ~contains(v,'unmy');
if any(sel), D.index(3) = find(sel,1); else D.index(3) = -1; end % myel axon diam
sel = contains(v,'diam') & contains(v,'ax') & contains(v,'unmy');
if any(sel), D.index(4) = find(sel,1); else D.index(4) = -1; end % unmy axon diam
sel = contains(v,'width');
if any(sel), D.index(5) = find(sel,1); else D.index(5) = -1; end % myelin width




return
function nerve = get_nerve_outline(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

% Gather data or parse inputs for spatial layout and types 
if any(named('-anat')), nerve = get_('-anat');
  if ischar(nerve) % filename
    nerve = mesh.insert_gmsh_fascicles('-info','-file',nerve,varargin{:});
  end
else nerve = mesh.insert_gmsh_fascicles('-info','-check-units');  
end
if isfield(nerve,'splines'), nerve = nerve.splines; end

% Get fasicle outlines
if iscell(nerve.outline), nerve.fascicles = nerve.outline{1}; 
else nerve.fascicles = nerve.outline; 
end

if isfield(nerve,'Attributes') && isfield(nerve,'Children') % XML leftovers 
    nerve = rmfield(nerve,{'Name','Attributes','Data','Children','Type'});
end

if isfield(nerve,'splines'), nerve = rmfield(nerve,'splines'); end

if max(abs(nerve.fascicles(:))) > 10 && ~any(named('-units-no-c'))
  nerve.coeffs = nerve.coeffs / 1e3;
  nerve.outline = nerve.outline / 1e3;
  nerve.fascicles = nerve.fascicles / 1e3;
  warning('ViNERS:make_axon_population:inputUnits', ...
          'Converting input units from µm to mm.')
end




%% CORE functions accomplish two aims: 
% (1) for 'marginal histogram' specifications, we need to compose a sample
%     population which generates the observed margainal histograms. 
function pop = generate_from_histograms(D, nom, types )

opts = optimoptions('fmincon','MaxFunctionEvaluations',1e5,'algorithm','sqp');

isa_C_type = contains(upper(types.Type),'UNMY');
A = sum(types.Count(~isa_C_type)); % Myelinated fibre count
C = sum(types.Count(isa_C_type)); % Unmyelinated fibre count
v_ = @(x) reshape(x,[],1);

pop = struct; 

% based on {nom}, this initialises:
% pop.fibre_diam, pop.fibre_gratio, pop.unmyelinated_diam

for f = 1:size(nom,1) % Initialise distributions

%       y = D.(nom{f,1}).y; 
%       lhs = diag(ones(numel(y)+1,1))/2 + diag(ones(numel(y),1),1)/2;
%       lhs(end,:) = []; 
%       y = pinv(lhs)*y
%%

    [y,vp] = unique(cumtrapz(D.(nom{f,1}).y));
    x = v_(D.(nom{f,1}).x(vp));      
    % x = max(0,[x;2*x(end)-x(end-1)]-median(diff(x))/2); % x = bin edges  
    % y = conv([0;y/max(y);1],[1;1]/2,'valid'); % [0; y/max(y)];
    % while y(2) == 0, y(1) = []; x(1) = []; end

    pop_ = @(xx,yy,nn) reshape(interp1(yy,xx,linspace(0,1,nn) ...
                    ,'linear','extrap'),[],1); 
    if strncmp(nom{f,1},'unmyelinated',3)
         pop.(nom{f,2}) = pop_(x,y,C);
    else pop.(nom{f,2}) = pop_(x,y,A);
    end

    if 0
      %% Debug visualisation

      clf, hold on %#ok<UNRCH>
      bar(D.(nom{f,1}).x,D.(nom{f,1}).y,0.95,'FaceColor',[.7 .7 .7],'EdgeColor','none')
      plot(D.(nom{f,1}).x, hist(pop.(nom{f,2}),D.(nom{f,1}).x)/numel(pop.(nom{f,2})),'-','LineWidth',1.2)

    end

    %% In case we don't do an optimiser later, do a 1D optimiser now
    n = numel(pop.(nom{f,2}));
    pop.(nom{f,2}) = smooth_exact_hist(D.(nom{f,1}), n);
end


if isfield(pop,'fibre_axon')
  if ~isfield(pop,'fibre_diam') && isfield(pop,'fibre_gratio')

    [~,shuffle] = sort(rand(size(pop.fibre_gratio)));
    pop.fibre_gratio = pop.fibre_gratio(shuffle);
    [~,shuffle] = sort(rand(size(pop.fibre_axon)));
    pop.fibre_axon = pop.fibre_axon(shuffle);

    mgt0 = @(x) min(x(x>0)); % QC check, remove negative values 
    pop.fibre_axon(pop.fibre_axon <= 0) = mgt0(pop.fibre_axon);
    pop.fibre_gratio(pop.fibre_gratio <= 0) = mgt0(pop.fibre_gratio);
    pop.fibre_diam =  pop.fibre_axon ./ pop.fibre_gratio;
  else
    y = (pop.fibre_axon > pop.fibre_diam); 
    d = pop.fibre_axon(y); 
    pop.fibre_axon(y) = pop.fibre_diam(y); 
    pop.fibre_diam(y) = d; 
  end
else        
  [~,shuffle] = sort(rand(size(pop.fibre_gratio)));
  pop.fibre_gratio = pop.fibre_gratio(shuffle);
  [~,shuffle] = sort(rand(size(pop.fibre_diam)));
  pop.fibre_diam = pop.fibre_diam(shuffle);

  mgt0 = @(x) min(x(x>0)); % QC check, remove negative values    
  pop.fibre_diam(pop.fibre_diam <= 0) = mgt0(pop.fibre_diam);
  pop.fibre_gratio(pop.fibre_gratio <= 0) = mgt0(pop.fibre_gratio);
  pop.fibre_axon = pop.fibre_gratio .* pop.fibre_diam;

end

%%
clear pdf_ pop_ f y vp nom  
if 0, mk_population_fitusing_ga, end %#ok<UNRCH>

lhs_AD = []; % LHS for axon diameter 
lhs_FD = []; % LHS for fibre diameter 
lhs_GR = []; % LHS for g ratio
% lhs_MW = []; % LHS for myelin width

graph_roi = []; % polyshape vertices for each ROI

id = fieldnames(D);
id(strcmp(id,'index')) = [];
var_ = @(x) D.(id{D.index(x)}); % g-ratio, fibre diam, axon diam, unmyel, width

if D.index(2) > 0, v2 = var_(2); else v2.x = 1; end
for fd = 1:numel(v2.x)
 if D.index(3) == -1 % missing axonDiameter
    warning('ViNERS:make_axon_population:noAxonDiameter', ...
            'missing axon (inner) diameter data, cannot estimate joint distribution of myelinated fibre sizes')
     break, 
 elseif D.index(2) == -1 % missing FibreDiameter
    warning('ViNERS:make_axon_population:noFibreDiameter', ...
            'missing fibre (outer) diameter data, cannot estimate joint distribution of myelinated fibre sizes')
    break, 
 end
 fd_r = v2.x(fd) + [-0.5 0.5] * mean(diff(v2.x));
 % if fd_r(1) == 0.5, fd_r(1) = 0; end % extend to zero

 v3 = var_(3); 

 for ad = 1:numel(v3.x)

  ad_r = v3.x(ad) + [-0.5 0.5] * median(diff(v3.x));
  % if ad == 1, ad_r(1) = 0; end % extend to zero
  if D.index(1) > 0
    gr_ok = [ad_r(1)/fd_r(2) min(1,ad_r(2)/fd_r(1))];
    v1 = var_(1);
    for gr = 1:numel(v1.x)

      gr_r =  v1.x(gr) + [-0.5 0.5]*median(diff(v1.x));
      if gr_r(1) < min(v1.x), gr_r(1) = 0; end
      if gr_r(2) > max(v1.x), gr_r(2) = 1; end
      if gr_r(1) - gr_ok(2) > -1e-8, continue, end
      if gr_r(2) - gr_ok(1) < 1e-8, continue, end

      graph_roi = [graph_roi; fd ad gr]; %#ok<AGROW>
      n = size(graph_roi,1);
      lhs_AD(ad,n) = 1; %#ok<AGROW>
      lhs_FD(fd,n) = 1; %#ok<AGROW>
      lhs_GR(gr,n) = 1; %#ok<AGROW>
    end
  elseif D.index(5) > 0
    mw_ok = [max(0,fd_r(1)-ad_r(2)) fd_r(2)-ad_r(1)];
    v5 = var_(5);
    for mw = 1:numel(v5.x)

      mw_r =  v5.x(mw) + [-0.5 0.5] * median(diff(v5.x));
      if mw_r(1) < min(v5.x), mw_r(1) = 0; end
      if mw_r(2) > max(v5.x), mw_r(2) = max(v5.x); end
      if mw_r(1) > mw_ok(2), continue, end
      if mw_r(2) < mw_ok(1), continue, end

      graph_roi = [graph_roi; fd ad mw]; %#ok<AGROW>
      n = size(graph_roi,1);
      lhs_AD(ad,n) = 1; %#ok<AGROW>
      lhs_FD(fd,n) = 1; %#ok<AGROW>
      lhs_GR(mw,n) = 1; %#ok<AGROW>
    end
  else
    if min(ad_r) > min(fd_r), continue, end
    graph_roi = [graph_roi; fd ad 0 1]; %#ok<AGROW>      
    n = size(graph_roi,1);
    lhs_AD(ad,n) = 1; %#ok<AGROW>
    lhs_FD(fd,n) = 1; %#ok<AGROW>
  end
 end
end

lhs = [lhs_FD; lhs_AD; lhs_GR]; 
if isempty(lhs), rhs = []; 
elseif D.index(1)>0 rhs = [v2.y; v3.y; v1.y]; % D.fibreDiam.y; D.axonDiam.y; D.g_ratio.y];
elseif D.index(5)>0 rhs = [v2.y; v3.y; v5.y]; % D.fibreDiam.y; D.axonDiam.y; D.sheathWidth.y];
else                rhs = [v2.y; v3.y]; % D.fibreDiam.y; D.axonDiam.y;
end

clear ad fd gr mw fd_r ad_r gr_r gr_ok mw_ok gr_rp 

if ~isempty(lhs) % 3 or more margainal histograms 

    est_init = pinv(lhs) * rhs;
    est_init(est_init < 0 ) = 0;
    est_init = est_init / sum(est_init(:));

    %%

    warning('off','MATLAB:nearlySingularMatrix')

    %                 @func, x0,           lhs*x = rhs, LB, UB, ... 
    est_fit = fmincon(@norm,est_init,[],[],lhs,rhs, ...
                      0*est_init,0*est_init+max(rhs),[],opts);
    warning('on','MATLAB:nearlySingularMatrix')
    est_fit(est_fit < 0) = 0; 


    infeasible = v2.x(graph_roi(:,1))< v3.x(graph_roi(:,2)); 
    est_fit(infeasible) = 0; 

    nr = numel(est_init);

    clf
    plot(1:nr,est_fit,1:nr,est_init,1:nr,pinv(lhs) * rhs)
    tools.tidyPlot
    rhs_out = lhs * est_fit; %#ok<NASGU>

    %% Generate populations from histograms 

    pop.fibre_diam = [];
    pop.fibre_axon = []; 
    pop.fibre_gratio = []; 

    cat = floor(est_fit*A);    
    cf = cumsum(mod(est_fit*A,1));
    c2 = arrayfun(@(x) find(cf>=x,1), rand(A-sum(cat),1)*cf(end));    
    vg = true(size(cat));
    for ii = 1:numel(c2), cat(c2(ii)) = cat(c2(ii)) + 1; end

    for pp = 1:numel(est_fit)
      if cat(pp) <= 0, continue, end

      % graph_roi = [graph_roi; fd ad gr]; %#ok<AGROW>
      fd_r = v2.x(graph_roi(pp,1)) + [-0.5 0.5] * median(diff(v2.x));
      % if fd_r(1) == 0.5, fd_r(1) = 0; end

      ad_r = v3.x(graph_roi(pp,2)) + [-0.5 0.5] * median(diff(v3.x));
      % if ad_r(1) == 0.5, ad_r(1) = 0; end

      xy = rand(10*A,2) .* [diff(fd_r) diff(ad_r)] + [fd_r(1) ad_r(1)];


      if D.index(1) > 0
        gr_r =  v1.x(graph_roi(pp,3)) + [-1 1] * median(diff(v1.x))/2;
        if gr_r(1) < min(v1.x), gr_r(1) = 0; end
        if gr_r(2) > max(v1.x), gr_r(2) = 1; end

        gr = xy(:,2) ./ xy(:,1); 
        ok = gr >= gr_r(1) & gr <= gr_r(2);
      elseif D.index(5) > 0
        mw_r =  v5.x(graph_roi(pp,3)) + [-1 1] * median(diff(v5.x))/2;
        if mw_r(1) < min(v5.x), mw_r(1) = 0; end
        if mw_r(2) > max(v5.x), mw_r(2) = max(v5.x); end
        mw = xy(:,1) - xy(:,2); 
        gr = xy(:,2) ./ xy(:,1);
        ok = mw >= mw_r(1) & mw <= mw_r(2);
      else gr_r = [0 1];
        gr = xy(:,2) ./ xy(:,1); 
        ok = gr >= gr_r(1) & gr <= gr_r(2);
      end

      if sum(ok) >= cat(pp)

          ok = find(ok);
          ok = ok(1:cat(pp));

          pop.fibre_diam = [pop.fibre_diam; xy(ok,1)];
          pop.fibre_axon = [pop.fibre_axon; xy(ok,2)];
          pop.fibre_gratio = [pop.fibre_gratio; gr(ok)];

      elseif any(ok) && D.index(1) >0 % tighten in bounds 

          fd_r = [min(xy(ok,1)) max(xy(ok,1))];
          ad_r = [min(xy(ok,2)) max(xy(ok,2))];
          xy = rand(A,2) .* [diff(fd_r) diff(ad_r)] + [fd_r(1) ad_r(1)];
          gr = xy(:,2) ./ xy(:,1); 
          ok = gr >= gr_r(1) & gr <= gr_r(2);

          if sum(ok) >= cat(pp)

            ok = find(ok);
            ok = ok(1:cat(pp));

            pop.fibre_diam = [pop.fibre_diam; xy(ok,1)];
            pop.fibre_axon = [pop.fibre_axon; xy(ok,2)];
            pop.fibre_gratio = [pop.fibre_gratio; gr(ok)];

          else error make_more_data
          end

      elseif D.index(1) < 0 % ~isfield(D,'g_ratio')

        if fd_r(1)-ad_r(2) < mw_r(1) % bottom right corner
          fd_r = [mw_r(1)+ad_r(1) fd_r(2)];
          ad_r = [ad_r(1) fd_r(2)-mw_r(1)];
        else
          fd_r = [fd_r(1) mw_r(2)+ad_r(2)];
          ad_r = [fd_r(1)-mw_r(2) ad_r(2)];
        end

        xy = rand(A,2) .* [diff(fd_r) diff(ad_r)] + [fd_r(1) ad_r(1)];
        mw = xy(:,1) - xy(:,2); 
        gr = xy(:,2) ./ xy(:,1);
        ok = mw >= mw_r(1) & mw <= mw_r(2);

        if sum(ok) >= cat(pp)

          ok = find(ok);
          ok = ok(1:cat(pp));

          pop.fibre_diam = [pop.fibre_diam; xy(ok,1)];
          pop.fibre_axon = [pop.fibre_axon; xy(ok,2)];
          pop.fibre_gratio = [pop.fibre_gratio; gr(ok)];

        else error make_more_data
        end

        % clf, hold on    
        % plot(fd_r([1 1 2 2 1]),ad_r([1 2 2 1 1]),'k-')
        % px = [fd_r(1) ad_r(2)+mw_r(2) nan ad_r(1)+mw_r(1) fd_r(2)];
        % py = [fd_r(1)-mw_r(2) ad_r(2) nan ad_r(1) fd_r(2)-mw_r(1)];
        % plot(px(:),py(:),'.r-')
        % plot(xlim,xlim,'-','Color',[0 0 0 0.3])
        % tools.tidyPlot, axis equal
        % xlabel('FD'),ylabel('AD')
      elseif all(gr > 1)
        error('%d axons requested with infeasible geometry',cat(pp))
      else
        % Trying a new polyshape approach
        
        if D.index(1) > 0
          gr_r =  v1.x(graph_roi(pp,3)) + [-1 1] * median(diff(v1.x))/2;
          if gr_r(1) < min(v1.x), gr_r(1) = 0; end
          if gr_r(2) > max(v1.x), gr_r(2) = 1; end
        else % if D.index(5) > 0
          vg(pp) = false; 
          continue
          % error todo_tiny_ROI_wo_G_ratio_hist
        end
        
        ps1 = polyshape(fd_r([1 1 2 2 1]), ad_r([1 2 2 1 1]));
        ps2 = polyshape(fd_r([1 1 2 2 1]), fd_r([1 1 2 2 1]) .* gr_r([1 2 2 1 1]));
        ps2 = intersect(ps1,ps2); 
        
        ps2.Vertices(end+1,:) = mean(ps2.Vertices);
        n2a = mod(0:cat(pp)-1,size(ps2.Vertices,1))+1;
        
        ps_gr = ps2.Vertices(n2a,2) ./ ps2.Vertices(n2a,1); 
        
        pop.fibre_diam = [pop.fibre_diam; ps2.Vertices(n2a,1)];
        pop.fibre_axon = [pop.fibre_axon; ps2.Vertices(n2a,2)];
        pop.fibre_gratio = [pop.fibre_gratio; ps_gr];
        
        % ps = intersect(ps, pr);
        % Try a new approach using polyshape

        warning('ViNERS:make_axon_population:tinyroi', ...
                'Injecting %d points into tiny G-ratio ROI',cat(pp))          
      end
    end

    if ~all(vg) % any "tiny ROI" issues to fix

      error todo_scatterBS

    end

else    
    pop.fibre_diam = v_(pop.fibre_diam);
    pop.fibre_axon = v_(pop.fibre_axon);
    pop.fibre_gratio = v_(pop.fibre_gratio);
    pop.unmyelinated_diam = v_(pop.unmyelinated_diam);
end


swap = (pop.fibre_axon > pop.fibre_diam);
pfa = pop.fibre_axon(swap);
pop.fibre_axon(swap) = pop.fibre_diam(swap);
pop.fibre_diam(swap) = pfa; 
  
  
% (1) This step is not necessary if we have (axon,fibre) diameter pairs
%     for the given nerve. Instead we sample the supplied population to
%     generate a new axon population
function pop = sample_from_population(D, pop, types )

isa_C_type = contains(upper(types.Type),'UNMY');
A = sum(types.Count(~isa_C_type)); % Myelinated fibre count
C = sum(types.Count(isa_C_type)); % Unmyelinated fibre count
named = evalin('caller','named'); 


if any(named('-keep'))
  fprintf('Keeping %d/%d myelinated axons\n', numel(pop.fibre_diam), A)
  A = numel(pop.fibre_diam); 
elseif numel(pop.fibre_diam) >= A

  [~,subset] = sort(rand(size(pop.fibre_diam))); 

  pop.fibre_diam = pop.fibre_diam(subset(1:A),:);
  pop.fibre_axon = pop.fibre_axon(subset(1:A),:);
  pop.fibre_gratio = pop.fibre_gratio(subset(1:A),:);    

elseif numel(pop.fibre_diam) < A
  nP = numel(pop.fibre_diam); 

  superset = ceil(nP*rand(1,A-nP));
  pop.fibre_diam = pop.fibre_diam([1:end superset],:);
  pop.fibre_axon = pop.fibre_axon([1:end superset],:);
  pop.fibre_gratio = pop.fibre_gratio([1:end superset],:);    
end
  
  if any(named('-keep'))
    fprintf('Keeping %d/%d unmyelinated axons\n', numel(pop.unmyelinated_diam), C)
    C = numel(pop.unmyelinated_diam); 
  elseif isempty(pop.unmyelinated_diam)
      
    if D.index(4) < 0 && C>0
      error('Missing unmyelinated axon diamter histogram')
    end
      
    id = fieldnames(D);
    id(strcmp(id,'index')) = [];
    var_ = @(x) D.(id{D.index(x)});      
    pop.unmyelinated_diam = smooth_exact_hist(var_(4), C);
      
  elseif numel(pop.unmyelinated_diam) >= C
    
    [~,subset] = sort(rand(size(pop.unmyelinated_diam))); 
    pop.unmyelinated_diam = pop.unmyelinated_diam(subset(1:C),:);
    pop.unmyelinated_xy = pop.unmyelinated_xy(subset(1:C),:);
    
  elseif numel(pop.unmyelinated_diam) < C
    nP = numel(pop.unmyelinated_diam); 
    superset = ceil(nP*rand(1,C-nP));
    pop.unmyelinated_diam = pop.unmyelinated_diam([1:end superset],:);
    pop.unmyelinated_xy = pop.unmyelinated_xy([1:end superset],:);
  end
  
  if any(named('-keep')) % adjust axon types to reflect Haveton tables 
    atc = types.Count(~isa_C_type);    
    types.Count(~isa_C_type) = round(atc / sum(atc) * A);
    
    atc = types.Count(isa_C_type);
    types.Count(isa_C_type) = round(atc / sum(atc) * C);
  end
  
return

% (2) the requested population needs to be laid out into one or more nerve
%     fascicles. 
function [pop,sam] = arrange_data_sample(pop, nG) % IN-CONTEXT function

named = evalin('caller','named'); 
get_ = evalin('caller','get_'); 

axon_types = evalin('caller','axon_types'); 
nerve = evalin('caller','nerve'); 
D = evalin('caller','D');

sam = struct;

if any(named('-remake')), sourceFile = get_('-remake'); 
  %%  
  sourceMap = load(sourceFile);
  
  if isfield(sourceMap,'axons') && isfield(sourceMap,'axon_subtype') % Curated axons file    
    sourceMap.pop = sourceMap.axons;   
    sourceMap.sam = sourceMap.axon_subtype;
  end
  
  gdata = [pop.fibre_diam pop.fibre_axon]; 
  
  if isfield(sourceMap,'sam')
    error TODO_from_SAM_field
  else
    %%
    sam.A_axon = []; 
    
    is_myel = [sourceMap.pop.myelinated];
    myel_st = cat(1,sourceMap.pop(is_myel).size_sample);
    myel_fd = cat(1,sourceMap.pop(is_myel).fibre_diam);
    myel_ad = cat(1,sourceMap.pop(is_myel).axon_diam);

    nG = max(myel_st);
    the = @(y,grp) arrayfun(@(g) median(y(grp == g)), 1:nG)'; 
    
    sam.A_diam = the(myel_fd,myel_st);
    sam.A_inner = the(myel_ad,myel_st);
    [~,sam.A_axon] = arrayfun(@(u) min((sam.A_diam-pop.fibre_diam(u)).^2 + ...
                                       (sam.A_inner-pop.fibre_axon(u)).^2 ), ...
                                           (1:numel(pop.fibre_axon))');
    
    unmy_st = cat(1,sourceMap.pop(~is_myel).size_sample);
    unmy_fd = cat(1,sourceMap.pop(~is_myel).fibre_diam);    
    
    
    sam.C_diam = the(unmy_fd,unmy_st);    
    sam.C_diam(isnan(sam.C_diam)) = []; 
    [~,sam.C_axon] = arrayfun(@(d) min(abs(sam.C_diam-d)), ...
                                           pop.unmyelinated_diam);                                         
  end
else
  %% Generate De Novo axon arrangements

  gdata = [pop.fibre_diam pop.fibre_axon]; 
  if any(named('-zscore')) % not recommended, doesn't converge as well with n
    gdata = zscore([pop.fibre_diam   pop.fibre_gratio]); 
  end

  opts = statset('MaxIter',1000);
  gm_fun = @(x,k)kmeans(x, k, 'replicate',5,'Options',opts);    

  sam.A_axon = gm_fun(gdata,nG);
  sam.A_diam = arrayfun(@(g) median(gdata(sam.A_axon == g,1)), 1:nG); 
  sam.A_inner = arrayfun(@(g) median(gdata(sam.A_axon == g,2)), 1:nG); 
  % sam.gratio = arrayfun(@(g) median(gdata(sam.A_axon == g,2)), 1:nG); 

  if ~isfield(pop,'unmyelinated_diam')

    warning('ViNERS:make_axon_population:missingCFibres', ...
            'no C-fibre data found, ignoring...')

    sam.C_diam = 0.25;
    sam.C_axon = [];
  else
    sam.C_diam = quantile(pop.unmyelinated_diam, 2*nG+1);
    sam.C_diam(1:2:end) = [];
    [~,sam.C_axon] = arrayfun(@(t) min(abs(t-sam.C_diam)), pop.unmyelinated_diam);
  end
  clear gdata gm_fun opts
end

%% Generate random axon arrangements 
nA = sum(axon_types.Count); 

nerve.f_bounds      = min(permute(min(nerve.fascicles,[],1), [2 3 1]),[],2)';
nerve.f_bounds(2,:) = max(permute(max(nerve.fascicles,[],1), [2 3 1]),[],2)';

% Pass 1 - random scatter to get baseline fill ratio
xy = rand(nA,2) .* [nerve.f_bounds(2) - nerve.f_bounds(1) ...
                    nerve.f_bounds(4) - nerve.f_bounds(3)] + ...
                   [nerve.f_bounds(1)   nerve.f_bounds(3)];

nF = size(nerve.fascicles,3);
f_id = arrayfun(@(n) in_loop(xy,nerve.fascicles(:,:,n)), 1:nF,'Unif',0);
f_id = [f_id{:}]; 

in_a_bv = false(size(f_id(:,1))); % In a blood vessel?
if iscell(nerve.outline)
  for bb = 1:size(nerve.outline{2},3)
    in_a_bv = in_a_bv | in_loop(xy,nerve.outline{2}(:,:,bb));
  end
end

% Pass 2 - overgenerate to get 1.1x as many as we should need 
n_ratio = mean(any(f_id,2) & ~in_a_bv); 
n_lloyd = 0; % I actually found that random points (no relaxation) better
             % matched real EM than relaxed points 
if any(named('-n-l')), n_lloyd = get_('-n-l'); end             


while true % might need multiple passes 
  xy = rand(ceil(1.1*nA/n_ratio),2) .* [nerve.f_bounds(2) - nerve.f_bounds(1) ...
                                        nerve.f_bounds(4) - nerve.f_bounds(3)] + ...
                                       [nerve.f_bounds(1)   nerve.f_bounds(3)];
  for lloyd_relax = 1:n_lloyd % make more uniform 
     [cp,vp] = voronoin(xy); 
      cp(1,:) = NaN;    
      xy = [cellfun(@(p) nanmedian(cp(p,1)), vp) ...
            cellfun(@(p) nanmedian(cp(p,2)), vp)];

      % h.XData = xy(:,1); h.YData = xy(:,2);
  end

  if n_lloyd == 0, [cp,vp] = voronoin(xy);   end
  
  f_id = arrayfun(@(n) in_loop(xy,nerve.fascicles(:,:,n)), 1:nF,'Unif',0);
  f_id = [f_id{:}];  ok = any(f_id,2);

  in_a_bv = false(size(f_id(:,1))); % repeat "in a blood vessel"
  if iscell(nerve.outline)
    for bb = 1:size(nerve.outline{2},3)
      in_a_bv = in_a_bv | in_loop(xy,nerve.outline{2}(:,:,bb));
    end
  end
  ok(in_a_bv) = 0; 

  vp(~ok,:) = []; xy(~ok,:) = []; f_id(~ok,:) = []; % Remove extrafascicular
  [f_id,~] = find(f_id'); % convert to #
  f_id = reshape(f_id,[],1); % column vector

  if size(xy,1) < nA % decrease n_ratio and try again
    if n_ratio < 0.01
      error('Axon location generation failed to produce axons contained in fascicles'), 
    end
    n_ratio = n_ratio * 0.8; continue
  else xy = xy(1:nA,:); f_id = f_id(1:nA,:); vp = vp(1:nA); break
  end
end

uatt = upper(axon_types.Type);

isa_C_type = contains(uatt,'UNMY');
isa_aff = contains(uatt,'AFF') | contains(uatt,'SENS');

nAA = sum(axon_types.Count(~isa_C_type & isa_aff));
nAE = sum(axon_types.Count(~isa_C_type & ~isa_aff));
nCA = sum(axon_types.Count(isa_C_type & isa_aff));
nCE = sum(axon_types.Count(isa_C_type & ~isa_aff));

perimeter_fcn = @(v) sum(sqrt(sum((cp(v,:)-cp(v([2:end 1]),:)).^2,2)));
[vcp,is_myelin] = sort(cellfun(perimeter_fcn, vp),'descend'); %#ok<ASGLU>
is_myelin = ismember((1:nA)',is_myelin(1:nAA+nAE));

is_afferent = is_myelin;
[~,aidx] = sort(rand(nAA+nAE,1));
is_afferent(is_myelin) = ismember(1:sum(is_myelin),aidx(1:nAA));

[~,aidx] = sort(rand(nCA+nCE,1));
is_afferent(~is_myelin) = ismember(1:sum(~is_myelin),aidx(1:nCA));

if isfield(D,'fibreEfferent')
  error apply_type_customisation_Adelta
end

if isfield(D,'unmyelinatedEfferent')
  %%
  [~,aidx] = sort(rand(size(pop.unmyelinated_diam))); 
  pop.unmyelinated_diam = pop.unmyelinated_diam(aidx); 
  
  unmy_aff = is_afferent(~is_myelin);
  
  x = mean(diff(D.unmyelinatedEfferent.x))/2;
  x = D.unmyelinatedEfferent.x-x; 
  x(end+1) = inf; x(1) = -inf;
  
  for ii = 1:numel(x)-1
    
    sel = (pop.unmyelinated_diam >= x(ii) & ...
           pop.unmyelinated_diam < x(ii+1));
   if ~any(sel), continue, end
   unmy_aff(sel) = linspace(0,1,sum(sel)) < D.unmyelinatedEfferent.y(ii); 
  end
  
  del = sum(unmy_aff) - axon_types(1,2);
  
  if del > 0
    aidx = find(unmy_aff);
    [~,sel] = sort(rand(size(aidx))); 
    unmy_aff(aidx(sel(1:del))) = 0;
  elseif del < 0
    aidx = find(~unmy_aff);
    [~,sel] = sort(rand(size(aidx))); 
    unmy_aff(aidx(sel(1:del))) = 1;
  end
  
  is_afferent(~is_myelin) = unmy_aff; 
    
 %%
end

pop.fibre_xy = xy(is_myelin,:);
pop.fibre_fascicle = f_id(is_myelin);
pop.fibre_afferent = is_afferent(is_myelin);

pop.unmyelinated_xy = xy(~is_myelin,:);
pop.unmyelinated_fascicle = f_id(~is_myelin);
pop.unmyelinated_afferent = is_afferent(~is_myelin);

% TODO : code to duplicate fascicle1 configuration for replicated fascicles
% (I wrote this , but somehow lost it possibly in a git pull or an experiment folder)

% pop = patch_fascIDs(pop,nerve); 

assignin('caller','xy',xy);
assignin('caller','f_id',f_id);
assignin('caller','is_myelin',is_myelin);
assignin('caller','is_afferent',is_afferent);

return

%%

clf %#ok<UNRCH>
scatter(pop.fibre_xy(:,1),pop.fibre_xy(:,2),10*pop.fibre_diam,pop.fibre_fascicle,'o')
colormap(lines(max(pop.fibre_fascicle))), axis equal, grid on
%%



%% Helper function 
function [D,nom] = parse_affeff_balance(D,nom)


warning('aff/eff not currently implemented. Also, this might crash')

nom(contains(nom(:,1),'aff'),2) = {'aff'};
nom(contains(nom(:,1),'eff'),2) = {'eff'};

return


%% are xy contained in the loop (alias inpolygon)
function is = in_loop(xy,loop)

is = inpolygon(xy(:,1),xy(:,2),loop(:,1),loop(:,2));
return

%%
is = false(size(xy(:,1)));
for ii = 1:size(xy,1)
    
    dx = sign(loop(:,1) - xy(ii,1));
    cx = find(dx ~= circshift(dx,[1 0])); 
    if isempty(cx), continue, end
    cn = 0*cx; 
        
    % clf, hold on, C = lines(7);
    % plot(path(:,1),path(:,2),'color',[0 0 0 0.3])
    % plot(xy(ii,1),xy(ii,2),'s','LineWidth',1.2,'Color',C(2,:))
    % plot(path(cx,1),path(cx,2),'.k','MarkerSize',8)    
    % plot([1 1]*xy(ii,1),ylim,'--','LineWidth',1.2,'Color',C(2,:))
    
    for pp = reshape(cx,1,[])        
        % plot(path(pp-[0 1],1),path(pp-[0 1],2),'-k','LineWidth',1.2)
        
        pp_01 = mod(pp-[1 2], size(loop,1))+1;
        
        if all(loop(pp_01,2) >= xy(ii,2)), cn(cx==pp)=1; continue
        elseif all(loop(pp_01,2) < xy(ii,2)),            continue
        end
        
        u = interp1(loop(pp_01,1),loop(pp_01,2),xy(ii,1));        
        cn(cx==pp) = (u > xy(ii,2));
    end
    
    is(ii) = mod(sum(cn),2)==1;
end

%% Alternate method of making a matching population (no joint dist data)
function out = smooth_exact_hist(dat, n)
    
x0 = reshape(dat.x,1,[]);
y0 = reshape(dat.y,1,[]);
dx = median(diff(dat.x)) / 2;
nB = round(n*y0); 

while sum(nB) < n % add points to match goal count    
    [~,s] = find(rand*sum(nB) >= cumsum(nB),1,'last');    
    if nB(s) == 0, [~,s] = max(nB); end
    if isempty(s), [~,s] = max(nB); end
    assert(nB(s) > 0)
    nB(s) = nB(s)+1; 
end

while sum(nB) > n % remove points to match goal count
    [~,s] = find(rand*sum(nB) >= cumsum(nB),1,'last');
    if nB(s) == 0, [~,s] = max(nB); end
    assert(nB(s) > 0)
    nB(s) = nB(s)-1; 
end

% mUse = (m0 + mFit)/2;

out = []; 

for ii = 1:numel(nB) % add points (linspace within bin) 
    if nB(ii) == 0, continue, end

    % u = rand(1,n); 
    % u = linspace(0,1,nB(ii)+1)*(mUse/y0); 
    % u = cumsum(u(1:end-1)) / sum(u(1:end-1)); 

    px = linspace(x0(ii)-dx,x0(ii)+dx,nB(ii)+1); 
    px = conv(px,[1 1]/2,'valid');
    out = [out; px'];
end

return

%% Visualise result and process

% px = [x0 - dx; x0 + dx];     
% py = [y0 - m0*dx; y0 + m0*dx];    
% clf, hold on    
% plot(px(:),py(:),'-','LineWidth',1.2)    
% py = [y0 - mFit*dx; y0 + mFit*dx];
% plot(px(:),py(:),'-','LineWidth',1.2)

clf, hold on, py = hist(out,dat.x); 
bar(dat.x,py,0.95,'FaceColor',[.7 .7 .7],'EdgeColor','none')    

py = rand(size(out)); 

for ii = 1:numel(nB)
sel = (out >= x0(ii)-dx) & (out < x0(ii)+dx);
py(sel) = py(sel) * (sum(sel)-0.1) + 0.05;
end

plot(out, py, '.w')
plot(dat.x, dat.y * n, 'LineWidth', 1.2,'Color','r')
    
%% Updated output code, supports more population diversity
function P = convert_output_structure(pop,sam)

types = pop.axon_table;

uatt = upper(types.Type);

isa_C_type = contains(uatt,'UNMY');
isa_aff = contains(uatt,'AFF') | contains(uatt,'SENS');

P = [];

for ty = 1:height(types) 
  
  this = struct;
  this.axon_type  = types.Type{ty};
  this.axon_model = types.Model{ty};
  this.afferent   = isa_aff(ty);
  this.myelinated = ~isa_C_type(ty);
  
  if isa_C_type(ty)
    
    sel = (pop.unmyelinated_afferent == isa_aff(ty));
    
    this.fibre_diam = pop.unmyelinated_diam(sel);
    this.axon_diam  = [];
    this.g_ratio    = [];
    
    this.axon_xy  = pop.unmyelinated_xy(sel,:);
    this.fascicle = pop.unmyelinated_fascicle(sel);
    this.size_sample = sam.C_axon(sel);
  else
    
    sel = (pop.fibre_afferent == isa_aff(ty));
        
    this.fibre_diam = pop.fibre_diam(sel);
    this.axon_diam  = pop.fibre_axon(sel);
    this.g_ratio    = pop.fibre_gratio(sel);
    
    this.axon_xy  = pop.fibre_xy(sel,:);
    this.fascicle = pop.fibre_fascicle(sel);
    this.size_sample = sam.A_axon(sel);
  end
  
  this.color = pop.axon_color(ty,:);

  if isempty(P), P = this; else P(end+1) = this; end
end

return
  
%% Save each population as a .csv or .xlsx table    
function export_tables(P,sam,nerve,ext,dest)

if ~isfolder(dest), mkdir(dest); end
fprintf('Saving tables (.%s) to %s\n',ext,tools.file('T',dest))
for ty = 1:numel(P)
  %% Organise into tables 
  T = rmfield(P(ty),{'axon_type','axon_model','color'});
  T.myelinated = repmat(T.myelinated,size(T.fascicle));
  T.afferent = repmat(T.afferent,size(T.fascicle));
  
  T.axon_x = T.axon_xy(:,1);
  T.axon_y = T.axon_xy(:,2);
  T.axon_xy = []; 
  
  v = fieldnames(T); 
  v(~cellfun(@(x) isempty(T.(x)),v)) = [];
  T = struct2table( rmfield(T,v) );
  
  T = T(:,[1:2 end-(3:-1:0) 3:end-4]);
      
  output_file = tools.file('get','%s/axons %s (%s) (%%d).%s','next',dest,P(ty).axon_type,P(ty).axon_model,ext);
  writetable(T,output_file)

end

%%
T = struct;
T.ID = [1:numel(sam.A_diam) 1:numel(sam.C_diam)]'; 
T.myelinated = [true(size(sam.A_diam)) false(size(sam.C_diam))]';
T.diameter_um = [sam.A_diam sam.C_diam]';
T.g_ratio  = [sam.gratio ones(size(sam.C_diam))]';
T.axon_count = [arrayfun(@(x) sum(sam.A_axon == x),1:numel(sam.A_diam)) ...
                arrayfun(@(x) sum(sam.C_axon == x),1:numel(sam.C_diam))]';

T = struct2table(T);              

output_file = tools.file('get','%s/axon sampling (%%d).%s','next',dest,ext);  
writetable(T,output_file)

if isempty(nerve), return, end

nerve.axon_color = cat(1,P.color);

output_file = tools.file('get','%s/../anatomy/anatomy (%%d).%s','next',dest,'mat');
if ~isfolder(fileparts(output_file)), mkdir(fileparts(output_file)); end
fprintf('Saving %s\n',tools.file('T',output_file))
save(output_file, '-struct','nerve');
  
%% Code for figure generation
function mk_hist_scatter_figure(D,pop)

id = fieldnames(D);
id(strcmp(id,'index')) = [];
var_ = @(x) D.(id{D.index(x)}); % g-ratio, fibre diam, axon diam, unmyel, width

r_ = @(x) reshape(x,1,[]);

%% Produce figure panels 
clf
set(gcf,'Position',get(gcf,'Position') .* [1 1 1.3 1])
xmax = round(quantile(pop.fibre_diam,0.998)); 

subplot(1,3,[1 2]), C = lines(7);
plot(pop.fibre_diam,pop.fibre_axon,'o','Color',C(2,:),'MarkerSize',4);
hold on,  axis equal, tools.tidyPlot, axis([-0.2 xmax -0.2 xmax])


if D.index(1) > 0
  v = var_(1); 
  x = xmax*r_(v.x);
  plot([0 xmax nan]'*(0*x+1),[0;1;nan]*x,'-','Color',[0 0 0 0.3])
elseif D.index(5) > 0  
  v = var_(5); 
  x = [0 ; r_(v.x) + 0.5 * mean(diff(v.x))];   
  px = [x 0*x+xmax nan*x]';
  py = [0*x xmax-x nan*x]';
  plot(px(:),py(:),'-','Color',[0 0 0 0.3])  
  
  v = var_(2); 
  px = [mean(v.x(1:2))*[1 1] nan xlim];
  py = [ylim nan mean(v.x(end-3:end-2))*[1 1]];
  plot(px,py,'k--','LineWidth',1.1)

end
x = 1:1:xmax;
plot([0 xmax nan]'*(0*x+1),[1;1;nan]*x,'-','Color',[0 0 0 0.3])

%%
bar_style = {'FaceColor',C(2,:),'EdgeColor',C(2,:), ...
             'LineWidth',1.1,'FaceAlpha',0.1};

subplot(3,3,3), cla, hold on % Fibre diameter
if D.index(2) > 0, v = var_(2);
  bar(v.x,v.y,1,'FaceColor',[.7 .7 .7],'EdgeColor','none')
      x = v.x;
else  x = 36;
end
[y,x] = hist(pop.fibre_diam,x); %#ok<*HIST>
bar(x,y/sum(y),1,bar_style{:})
tools.tidyPlot, xlabel('fibre diameter (µm)')

subplot(3,3,6), hold on

if D.index(3) > 0, v = var_(3);
 bar(v.x,v.y,1,'FaceColor',[.7 .7 .7],'EdgeColor','none')    
 x = v.x;
elseif D.index(2) > 0, v = var_(2); x = v.x;
else x = 36; 
end

[y,x] = hist(pop.fibre_axon,x); 
bar(x,y/sum(y),1,bar_style{:})     

tools.tidyPlot, xlabel('axon diameter (µm)')

subplot(3,3,9), hold on
if D.index(1) > 0, v = var_(1); 
  if isfield(v,'y')
    bar(v.x,v.y,1,'FaceColor',[.7 .7 .7],'EdgeColor','none')      
  end
  
  [y,x] = hist(pop.fibre_gratio,v.x); 
  bar(x,y/sum(y),1,bar_style{:})
  tools.tidyPlot, xlabel('g-ratio')
elseif numel(D.index) > 4 && D.index(5) > 0, v = var_(5);
    
  bar(v.x,v.y,1,'FaceColor',[.7 .7 .7],'EdgeColor','none')    
  y = hist(pop.fibre_diam - pop.fibre_axon,v.x); 
  bar(v.x,y/sum(y),1,bar_style{:})
  tools.tidyPlot, xlabel('sheath width (µm)')
end

return


%% Manually try to improve the histogram to population conversion (unused)
function manual_improve_hist_fit %#ok<DEFNU>

[~,idx] = sort(rand(A,1) - 0.5*linspace(0,1,A)'); 
axonDiam = pop.fibreDiam .* pop.g_ratio(idx); 

clf, C = lines(7);
h = plot(pop.fibreDiam,axonDiam,'o','Color',[.4 .4 .4],'MarkerSize',4);
hold on,  axis equal, tools.tidyPlot, axis([-0.2 8 -0.2 8])

x = 8*D.g_ratio.x;
plot([0 8 nan]'*(0*x+1),[0;1;nan]*x,'-','Color',[0 0 0 0.3])
x = 1:1:7;
plot([0 8 nan]'*(0*x+1),[1;1;nan]*x,'-','Color',[0 0 0 0.3])

h(2) = plot(0,0,'s','Color',C(2,:),'MarkerSize',10,'LineWidth',1.1);

axes('Position',[.2 .6 .3 .25])
bar(D.axonDiam.x,D.axonDiam.y,1,'FaceColor',[.4 .4 .4],'EdgeColor','none')
tools.tidyPlot, hold on

y = hist(axonDiam,D.axonDiam.x-0.5); %%#ok<HIST>
h(3) = plot(D.axonDiam.x,y/sum(y),'-+','Color',C(2,:),'LineWidth',1.1,'MarkerSize',10); 

swap = []; 
%%

while(true)
    
    [x,y,b] = ginput(1);
    [~,sel] = min((pop.fibreDiam-x).^2 + (axonDiam-y).^2);
    
    if isempty(b), break, end
    if b > 5, break, end
    if isempty(swap) || b == 3
        swap = sel;
        h(2).XData = h(1).XData(sel);
        h(2).YData = h(1).YData(sel);
    
    elseif ~isempty(swap)
        
        idx([swap sel]) = idx([sel swap]);
        axonDiam = pop.fibreDiam .* pop.g_ratio(idx); 
        h(1).YData = axonDiam;
        h(2).XData = [h(2).XData h(1).XData(sel)];
        h(2).YData = axonDiam([swap sel]);                
        
        y = hist(axonDiam,D.axonDiam.x-0.5);
        h(3).YData = y./sum(y);
        swap = []; 
    end
end
clear h





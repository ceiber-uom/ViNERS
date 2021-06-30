

function view_thresholds(varargin) % Adapted from EMBC fig. 4

named = @(v) strncmpi(v,varargin,length(v));
get_ = @(v) varargin{find(named(v))+1};
f_ = @(x) [x.folder filesep x.name];

% Choose simulation to display
if nargin > 0 && exist(varargin{1},'dir')
  stim_folder = varargin{1}; 
elseif nargin > 0 && strcmp(varargin{1},'?')
  stim_folder = uigetdir(tools.file('thresholds~/'));
else
  stim_folder = tools.file('thresholds~/');
  list = dir([stim_folder '*']);  
  list(cellfun(@(n) all(n=='.'),{list.name})) = [];
  if any([list.isdir]), list = list([list.isdir]); end
    
  if nargin > 1, sel = strcmp(varargin{1},{list.name});
  else           sel = false;
  end
  
  if any(sel), stim_folder = f_(list(find(sel,1)));
  elseif ~any([list.isdir]),     
    stim_folder = uigetdir(tools.file('stim~/')); % try here?
  elseif numel(list) == 1, stim_folder = f_(list);
  else    stim_folder = uigetdir(list(1).folder);
  end
end

%%

load_vars_ = @(s) load(s,'diameter','threshold','axon_index');
has_vars_ = @(s,varargin) ~isempty(whos('-file',s,varargin{:}));

if any(named('-ax')), ax_file = get_('-ax'); 
else ax_file = tools.file('get','axons~/axon*.mat','newest');
end


if ~any(named('-q'))
  
  fprintf('Viewing %s\nLoading %s\n',tools.file('T',stim_folder), ...
                                     tools.file('T',ax_file))
  
end

load(ax_file,'nerve','pop'); 

if ~exist('nerve','var')
  fasc_file = tools.INPUT_file('splines.dat'); 
  nerve = mesh.read_dat_file(fasc_file);
  nerve.fascicles = nerve.outline{1};
end
% pop.axon_color = pop.axon_color([3 4 1 2],:);
% nerve.fascicle = nerve.outline{1};

%% Set up fascicles to the left

figure(1), clf, G = @(v) [v v v]/10; 
outline_style = {'EdgeColor',G(2),'LineWidth',1.1};

nF = size(nerve.fascicles,3); 

if max(abs(nerve.fascicles(:))) > 10 && ~any(named('-units-um'))
  disp('CONVERTING FASCICLES TO MM')  
  nerve.fascicles = nerve.fascicles / 1e3;
  
end

for aa = 1:2
  subplot(3,2,2*aa-1), cla, hold on
  % old epineurium
  % fill(nerve.outline{3}(:,1), nerve.outline{3}(:,2), G(8.5), outline_style{:})
  
  bg = G(8.5); 
  
  for ff = 1:nF
    fill(nerve.fascicles(:,1,ff),nerve.fascicles(:,2,ff),bg,outline_style{:})
  end
  axis tight equal, tools.tidyPlot

end

m = tools.magma(110);
colormap(gcf,flipud(m(1:100,:)))


%% Scatterplots

v_ = @(x) reshape(x,[],1);  %#ok<NASGU>
I_stim = cell(4,1);

list = dir([stim_folder '/*.mat']);


axon_colors = []; 
axon_names = {}; 

for ff = 1:numel(list)
  
  f_id = str2double(regexp(list(ff).name,'(?<=fasc[^\d]+)\d+','match'));  
  stim = load_vars_(f_(list(ff)));
  
  if has_vars_(f_(list(ff)),'pop'), load(f_(list(ff)),'pop'); end

  is_my = pop.myelinated; 
  diam = stim.diameter;  
  if size(diam,1) ~= size(stim.threshold,1), 
      diam = diam(pop.fascicle == f_id,:);
  end
  xy = pop.axon_xy(pop.fascicle == f_id,:);
    
    
  subplot(3,2,3-2*is_my); hold on
  scatter(xy(:,1),xy(:,2),diam*3,stim.threshold,'o','LineWidth',1.1)  
    
  subplot(3,2,4-2*is_my); hold on
  scatter(diam,stim.threshold,diam*3,stim.threshold,'o','LineWidth',1.1)

  xlabel('fibre diameter (µm)')
  ylabel('threshold (µA)')
  
  if numel(pop.population_id) == numel(pop.fascicle)
    for ty = unique(pop.population_id)'      
      I_stim{ty} = [I_stim{ty}; stim.threshold(pop.population_id(pop.fascicle == f_id) == ty)];
    end, ty = unique(pop.population_id);
  else ty = pop.population_id;
      I_stim{ty(1)} = [I_stim{ty(1)}; stim.threshold];
  end
    
  axon_colors(ty,:) = pop.color;
  if iscell(pop.axon_model), axon_names(ty) = pop.axon_model;
  else                       axon_names(ty) = {pop.axon_model};
  end
end

most_of = 0.01;
if any(named('-q')), most_of = get_('-q'); end


fudge = cat(1,I_stim{:});
fudge = std(fudge(isfinite(fudge))) / 1e6; 

most_of = @(n) quantile(cat(1,I_stim{n}),[most_of 1-most_of]) + [-1 1]*fudge;

if ~isempty(I_stim{1}) || ~isempty(I_stim{2}) , cax = most_of([1 2]); 
 if all(isfinite(cax))
  subplot(3,2,1), caxis(cax)
  subplot(3,2,2), caxis(cax)
 end, tools.tidyPlot
end
if ~isempty(I_stim{3}), cax = most_of([3 4]); 
 if all(isfinite(cax))
  subplot(3,2,3), caxis(cax)
  subplot(3,2,4), caxis(cax)
 end, tools.tidyPlot
end

%%

I_stim{1}(I_stim{1} == 0) = NaN;

mk_sumplot = @(n) plot(sort(I_stim{n}),linspace(0,100,numel(I_stim{n})),'.', ...
                           'Color',axon_colors(n,:),'LineWidth',1.2);

mk_mnln = @(n) plot([1 1]*nanmedian(I_stim{n}),[0 50], 'Color', ...
                             [axon_colors(n,:) 0.4],'LineWidth',1.2);

subplot(3,3,7), cla reset, hold on, 
mk_sumplot(1); mk_sumplot(2); mk_mnln(1); mk_mnln(2);
tools.tidyPlot, % xlim([0 200])
xlabel('threshold (µA)')

subplot(3,3,[8 9]),cla reset,  hold on, 
mk_sumplot(1); mk_sumplot(2); mk_mnln(1); mk_mnln(2);
mk_sumplot(3); mk_sumplot(4); mk_mnln(3); mk_mnln(4);
tools.tidyPlot,%  xlim([200 1600])
xlabel('threshold (µA)')

%%

report_ = @(n) fprintf('%s: %0.1f, %0.1f = 50/95%%\n', ... 
              axon_names{n(1)},quantile(cat(1,I_stim{n}),[.5 .95])); 

report_(1)
report_(2)
report_([3 4])

% max(I_stim{2})
% mean(cat(1,I_stim{[3 4]}) < max(cat(1,I_stim{[1 2]})))


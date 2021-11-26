
function preview_spike(filename,varargin)
% preview the spike response to extracellular stimuluation at various 
% fixed levels as generated by models.nerve_stimulation
% (extracted from models.axon_population)
% 
% Options:
% 
% -xy [x y] : view axon closest to specified xy point
% -id [id]  : view axon ID
% (default) : view axon with median threshold
% -no-anat  : do not generate anatomy panel inset
% -no-elec  : hide electrodes in generated plot
% -stim [f] : use specified eidors_file for electrode layout
% -fasc [a] : use specified fascicle anatomy
% -sp       : generate one row per stimulus (subplots)
%             -spp [percent overlap] and -spc [row color] fine-tune this
% -cv       : show line indicating conduction velocity (to check units).
% 
% V0.2 CDE 18-Nov-2021

if nargin > 0 && filename(1) == '-'
    varargin = [{filename} varargin]; filename = ''; 
end

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if nargin == 0 || isempty(filename) || strcmp(filename,'?')
     filename = tools.file('get','stim~/s*');
     if isempty(filename), filename = tools.file('stim~/'); end
    [filename,fp] = uigetfile('*.mat','',filename);
     filename = [fp filename];
elseif any(filename == '~'), filename = tools.file(filename);
end

fprintf('reading %s\n', tools.file('T',filename))
d = load(filename); 

fibre_xy = d.pop.axon_xy(d.axon_index,:);

if any(named('-xy'))
  if any(named('-xyc')), xyc = mean(fibre_xy); 
  else xyc = get_('-xy');
  end
  [~,sel] = min(sum((fibre_xy-xyc).^2,2));
elseif any(named('-i')), sel = get_('-i'); 
else  
  [~,sel] = min(abs(d.threshold - median(d.threshold)));
end

%% Get other inputs wanted


if ~any(named('-no-e'))
  %%
  if any(named('-stim')), eidors_file = get_('-stim');
    if any(eidors_file == '~'), eidors_file = tools.file(eidors_file); end
  else eidors_file = tools.file('get','eidors~/stim*.mat');
  end

  fprintf('reading %s\n', tools.file('T',eidors_file))

  EM = load(eidors_file); 
  nE = size(EM.v_extracellular,2); 

  elec_y = [cellfun(@(e) min(EM.model.nodes(e,1)), {EM.model.electrode.nodes}); ...
            cellfun(@(e) max(EM.model.nodes(e,1)), {EM.model.electrode.nodes}) ];

else nE = 0; 
end



%%
if ~any(named('-no-a')) && ~any(named('-clf')), clf, end
cla reset, hold on
C = lines(7); 

do_subplots = any(named('-sp'));
if do_subplots
  do_sp_olap = 0.02; 
  if any(named('-spp')), do_sp_olap = get_('-spp')/100; end
  light_gray = [.9 .9 .9];
  if any(named('-spc')), light_gray = get_('-spc'); end
end

if ~isfield(d.stimulus,'current'), 
  d.stimulus.current = mean(d.threshold(isfinite(d.threshold)));
end

if ~any(size(d.stimulus.current) == 1)
  
  if any(named('-cid')), cid = get_('-cid');
  else warning('ViNERS:previewSpike:cid', ...
               'd.stimulus.current has multiple levels per stimulus, using row #1. use -cid to access other rows.')
      cid = 1;
  end
  d.stimulus.current = d.stimulus.current(cid,:); 
end

iMax = max(d.stimulus.current(:));
nStim = size(d.spike_counts,2);
if isfield(d.results(sel),'length')
    len_ = @(n) d.results(sel).length(n+1);
else
    % per-axon simulated length found here apparently
    length_to_mm = d.results(sel).summary(end);
    if isfield(d.results(sel),'n_nodes')
         len = d.results(sel).n_nodes;
    else len = max(d.results(sel).spikes.node{end});
    end
    len_ = @(n) ((n/len)-0.5 ) * length_to_mm  ;
end

iThreshold = d.threshold(sel); 
if ~isfinite(log10(iThreshold)), iThreshold = min(d.stimulus.current); end

color = interp1(log10([1; iThreshold; iMax]), ... 
                      [0 0 0; .1 .1 .1; C(7,:)], log10(d.stimulus.current)); 
for i_stim = 1:nStim
  
  if do_subplots
    subplot(nStim,1,i_stim), cla, hold on
  end
  
  t = d.results(sel).spikes.time{i_stim};
  x = d.results(sel).spikes.node{i_stim};  
  x = len_(x);
    
  plot(t,x,'.','Color',color(i_stim,:))
  
  [~,n] = min(x); 
  
  text(t(n), x(n)-mean(diff(x)), sprintf('%0.1f �A', ...
                        d.stimulus.current(i_stim)), ...
       'Horiz','center','vert','top','FontSize',7,'Color',color(i_stim,:))  
  
  if i_stim == 1
    
    px = d.stimulus.timing(max(abs(d.stimulus.pattern),[],2) > 1e-9);
    plot([0 range(px)],[1 1]*max(ylim),'k-')
  end
  
  if do_subplots
    if i_stim < nStim, set(gca,'XColor','none','YColor','none'), end    
    if mod(i_stim,2) == 0 && any(light_gray>0), set(gca,'Color',light_gray), end
    for ee = 1:nE
      plot(-0.02*[1 1], elec_y(:,ee),'-','Color',[.4 .4 .4],'LineWidth',1.5)
    % text(0,mean(elec_y(:,ee)),sprintf('E%d',ee),'Color',[.4 .4 .4],'horiz','left')
    end
  end

end

% plot([0 range(d.stimulus.timing)],[1 1]*max(ylim),'k-')
tools.tidyPlot

if do_subplots
  
  h = get(gcf,'children');
  
  set(h,'XLim',[-0.1 max([h.XLim])])
  set(h,'YLim',[min([h.YLim]) max([h.YLim])])
  for ee = 1:nE
      % plot(-0.02*[1 1], elec_y(:,ee),'-','Color',[.4 .4 .4],'LineWidth',1.5)
      text(-0.2,mean(elec_y(:,ee)),sprintf('E%d',ee),'Color',[.4 .4 .4],'horiz','right','Parent',h(end))
  end
  
  for ii = 1:numel(h)
    set(h(ii),'Position',get(h(ii),'Position')-[0 1 0 -2]*do_sp_olap)
  end
  set(gcf,'Children',flipud(get(gcf,'Children')))
  
  return
end

if any(named('-cv'))
  plot(d.results(sel).dt_dx*(1:len), len_(1:len),'--','Color',[C(3,:) 0.5])
end

%% Add annotations (not relevent to do_subplots)

if ~any(named('-no-e'))  
  %%
  for ee = 1:nE

    plot(-0.02*[1 1], elec_y(:,ee),'-','Color',[.4 .4 .4],'LineWidth',1.5)
    text(0,mean(elec_y(:,ee)),sprintf('E%d',ee),'Color',[.4 .4 .4],'horiz','left')
  end
end

if ~any(named('-no-a'))
  %%  
  h = legend('example','b','c','d','e','location','best'); 
  p = h.Position; delete(h);  
  axes('Position',p); 

  if any(named('-fasc')), nerve = get_('-fasc'); 
  else   
    axon_file = tools.file('sub~\axons\');
    if any(named('-root')), axon_file = get_('-root'); end
    axon_file = tools.file('get',[axon_file '/axons*.mat'],'newest');
    if isempty(axon_file)
      warning('axons.mat missing, please specify using -root [axon-file]')
      [fn,fp] = uigetfile('axons*.mat',[],tools.file('axons~\'));
      axon_file = [fp fn]; clear fp fn
    end
    fprintf('loading %s\n',tools.file('T',axon_file)) 
    load(axon_file,'nerve');
  end
  
  if range(nerve.fascicles(:)) > 10, nerve.fascicles = nerve.fascicles / 1e3; end
  if any(named('-dxy')), nerve.fascicles = nerve.fascicles + get_('-dxy'); end
  
  f_id = d.pop.fascicle; f_id = f_id(d.axon_index(sel)); 
  
  hold on
  plot(fibre_xy(:,1), fibre_xy(:,2),'.','Color',[.8 .8 .8])  
  plot(nerve.fascicles(:,1,f_id), nerve.fascicles(:,2,f_id), 'Color',[.3 .3 .3],'LineWidth',1.2)
  
  plot(fibre_xy(sel,1), fibre_xy(sel,2),'s','Color',C(7,:))  
  axis image, tools.tidyPlot
end
  
  



function preview_axons(filename,varargin)
% preview the axon arrangement in a group of fascicles 
% (extracted from models.axon_population)
% 
% -down [x]  : apply downsample factor for display
% -fig1 [f]  : plot anatomy in specified figure (Default: figure #1)
% -no-labels : do not write fascicle labels on figure
% -elec [info] : show electrodes (info is output of models.electrode_array)
% -no-fig2   : do not generate second figure 
% -fig2 [f]  : plot histograms in specified figure (Default: figure #2)
% -sampling  : show sample assignment (voronoi diagram)
% -no-patch  : do not fill in patch areas (voronoi diagram)
% -g-ratio   : plot g-ratio scatterplot instead of 
% 
% V0.3 CDE 18-Nov-2021


named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if nargin == 0, filename = tools.file('get','axons~/ax*.mat'); end
if contains(filename,'~'), filename = tools.file(filename); end
if ~exist(filename,'file')
  filename = tools.file('axons~/axons*.mat','-prompt');  
end

load(filename,'nerve','pop');

if ~exist('nerve','var')
  nF = max(pop.unmyelinated_fascicle);
  nerve.coeffs = zeros(2,0,nF);
  nerve.fascicles = zeros(0,2,nF);
  nerve.outline = zeros(0,2,nF);
end

if ~isfield(nerve,'fascicles')  
  nerve.fascicles = nerve.outline;
  assert(~iscell(nerve.outline))
end

nF = size(nerve.coeffs,3);
ds_factor = 1; % no downsampling

plot_g_ratio = 0; 
if any(named('-g-a')), plot_g_ratio = 2;
elseif any(named('-g')), plot_g_ratio = 1;
end
    

if any(named('-down')), ds_factor = get_('-down'); end

%%
if any(named('-fig')), figure(get_('-fig')); 
else figure(1), p = get(gcf,'Position'); 
  if all(p(3:4) == [560 420]), set(gcf,'Position',p .* [0.8 1 1.3 1]), end
end

clf, hold on, % C = lines(nF);

labels = {'Fascicle outline','unmyelinated afferent','unmylinated efferent', ...
            'myelinated afferent','myelinated efferent'};

if max(nerve.fascicles(:)) > 20 % mm to um
  nerve.fascicles = nerve.fascicles / 1000; 
end         
if iscell(nerve.outline)
  fill(nerve.outline{end}(:,1),nerve.outline{end}(:,2),[.9 .9 .9],'EdgeColor','none')
  labels = [{''} labels];
end

for ff = 1:nF
    fill(nerve.fascicles(:,1,ff),nerve.fascicles(:,2,ff),'w','EdgeColor',[.3 .3 .3],'LineWidth',1.2)
    
    if ~any(named('-no-l')) 
      text(mean(nerve.fascicles(:,1,ff)),mean(nerve.fascicles(:,2,ff)),sprintf('\\bfF%d',ff), ...
          'Color',[.7 .7 .7],'FontSize',26,'Horiz','center')
    end

    for pp = numel(pop):-1:1
    
        f_id = (pop(pp).fascicle == ff); 
        xy = pop(pp).axon_xy;
    
        if pop(pp).myelinated
            style = {'o','MarkerSize',5, ...
                         'MarkerFaceColor',(pop(pp).color+1.5)/2.5', ...
                         'LineWidth',1.2};
        else style = {'.','MarkerSize',10};
        end
        
        f_id(f_id) = mod(1:sum(f_id), ds_factor) == 0;
        plot(xy(f_id,1),xy(f_id,2),style{:},'color',pop(pp).color)
    end
end

if any(named('-e')), e = get_('-e'); 
  if isfield(e,'info'), e = e.info; end
  
  y2x = range(ylim) / range(e.ElectrodePositions(:,1));
  
  for cc = 1:size(e.ElectrodePositions,1)
    
    xy = [e.ElectrodePositions(cc,[3 1]) 0 0] + [-1 -1 2 2]/2 .* ...
          e.ElectrodeDimensions(e.ElectrodeTypeIndex(cc),[3 1 3 1]);
    
    xy(2) = xy(2) - max(e.ElectrodePositions(:,1));
    xy([2 4]) = xy([2 4]) * y2x;
    
    rectangle('Position',xy,'FaceColor',[.6 .6 .6],'EdgeColor',[.4 .4 .4])
    
  end
end

axis image, tools.tidyPlot
title(regexprep(tools.file('T',filename),'([\\_\^])','\\$1'))
legend(labels{:},'location','best')

if any(named('-no-f')), return, end

%%

if any(named('-fig2')), figure(get_('-fig2')); 
else figure(2), p = get(gcf,'Position'); 
  if all(p(3:4) == [560 420]), set(gcf,'Position',[220 120 730 700]), end
end

%% Produce figure panels 
clf
xmax = ceil(quantile(cat(1,pop.fibre_diam),0.998)); 
bar_style = {'EdgeColor','none','FaceAlpha',0.5,'FaceColor'};

[~,x] = hist(cat(1,pop.fibre_diam),36); %#ok<HIST>

G = @(v) [v v v]/10; C = lines(7);

cla
for pp = find([pop.myelinated])  
    subplot(3,3,[4 8]), hold on,
    
    X = pop(pp).fibre_diam;
    Y = pop(pp).axon_diam; 
    
    if plot_g_ratio > 0, Y = Y./X; end
    if plot_g_ratio > 1, X = pop(pp).axon_diam; end
    
    plot(X,Y,'o','Color',pop(pp).color,'MarkerSize',4)     
  
end

if any(named('-sam'))
    %%
    
    xmax = ceil(quantile(cat(1,pop.fibre_diam),1)); 

    ax_l = axis; cla, hold on
    
    sel = [pop.myelinated]; 
    gid = cat(1,pop(sel).size_sample);
    
    [xy,id] = unique([cat(1,pop(sel).fibre_diam) ...
                      cat(1,pop(sel).axon_diam)],'rows');
    gid = gid(id);
    
    g_color = zeros(max(gid),4);
    
    [vert,cb] = voronoin(xy);
    neighbour = delaunayn(xy); 
    
    for ii = 1:size(cb,1)
        
      if ~any(g_color(gid(ii),:))
        
        these = find(gid == gid(ii));
        sel = any(ismember(neighbour,these),2);
        sel = setdiff(neighbour(sel,:), these); 
        sel = unique(gid(sel)); 
        
        cid = setdiff(1:7,g_color(sel,4)); 
        cid = cid(ceil(rand*numel(cid)));
        
        g_color(gid(ii),4) = cid(1);
        g_color(gid(ii),:) = [C(cid(1),:) cid(1)];
      end
      
      if any(named('-no-p')), continue, end
      
      [vi,~,vb] = unique(cb{ii}(cb{ii} ~= 1)); 
      vx = vert(vi,:); 
      
      if plot_g_ratio > 1,
              vx = fliplr(vx); vx(:,2) = vx(:,1) ./ vx(:,2); 
      elseif plot_g_ratio > 0, vx(:,2) = vx(:,2) ./ vx(:,1); 
      end
      
      patch('Faces',vb','Vertices',vx,'FaceColor', ...
                    g_color(gid(ii,:),1:3),'EdgeColor','none','FaceAlpha',0.2)
    end
    axis(ax_l)

    if plot_g_ratio > 1,
            xy = fliplr(xy); xy(:,2) = xy(:,1) ./ xy(:,2); 
    elseif plot_g_ratio > 0, xy(:,2) = xy(:,2) ./ xy(:,1); 
    end

    
    scatter(xy(:,1),xy(:,2),[],g_color(gid,1:3),'.')
    
    for gg = 1:max(gid)
        
        vi = [cb{gid == gg}]; 
        vert(vert > max(xlim)) = max(xlim);
        vx = mean(vert(vi,:));        
        ap = median(xy(gid == gg,:),1);
        ap(2) = ap(2)./ap(1);
        text(vx(1),vx(2),sprintf('#%d\n^{d=%0.2f,g=%0.2f}',gg,ap),'FontSize',8,'Color',g_color(gg,1:3))        
    end
    
    clear ax_l vi vb vx ap cid sel these 
end

%%

is_myel = find([pop.myelinated]);   
    
if any(is_myel)
  
    subplot(3,3,[4 8]),
    if plot_g_ratio == 0, axis equal, end
    tools.tidyPlot, axis([-0.2 xmax -0.2 xmax]), grid on
      
    h = flipud(get(gca,'children'));
    
    if plot_g_ratio == 0, % diagonal grid lines for g-ratio 
      gid = 0.1:0.1:1;
      plot([0;xmax;nan]*(0*gid+1),[0;xmax;nan]*gid,'-','Color',[0 0 0 0.3])
      text([xmax 0]*[1.01; -0.01],0.1*xmax,'g=0.1','Color',G(7))
      text([xmax 0]*[1.01; -0.01],0.9*xmax,'g=0.9','Color',G(7))
    end
    
    if ~any(named('-sam'))
      legend(h,pop(is_myel).axon_type,'location','nw','box','off')
    end

    xlabel('fibre diameter (µm)'),
    ylabel('axon diameter (µm)'),
    if plot_g_ratio > 0, ylim([-0.02 1]), ylabel('g ratio'), end
    if plot_g_ratio > 1, xlabel('axon diameter (µm)'), end
    
    subplot(3,3,[1 2]), hold on        
    if plot_g_ratio > 1,       
         y = hist(cat(1,pop(is_myel).axon_diam),x); %#ok<HIST>
    else y = hist(cat(1,pop(is_myel).fibre_diam),x); %#ok<HIST>
    end
    
    bar(x,y,1,bar_style{:},pop(is_myel(end)).color)
    
    xlim([0 xmax]), tools.tidyPlot, ylabel('myelinated')
    set(gca,'YColor',pop(is_myel(end)).color)
    ylim([0 max(ylim)])

    if any(~[pop.myelinated]), axes('position',get(gca,'Position')), end
else cla
end

is_unmy = setdiff(1:numel(pop),is_myel); 

if any(is_unmy)
  y = hist(cat(1,pop(is_unmy).fibre_diam),x); %#ok<HIST>
  bar(x,y,1,bar_style{:},pop(is_unmy(1)).color)

  xlim([0 xmax]), tools.tidyPlot,  ylabel('unmyelinated')
  set(gca,'YColor',pop(is_unmy(1)).color,'YAxisLocation','right')
  ylim([0 max(ylim)])
end

if any(is_myel), subplot(3,3,[6 9]),

    [y,x] = hist(cat(1,pop(is_myel).g_ratio),0:0.05:1); %#ok<HIST>
    barh(x,y,1,bar_style{:},pop(is_myel(end)).color)
    tools.tidyPlot, ylabel('g-ratio')
    set(gca,'YTick',0:0.1:1)
    xlabel('axon count')


end

%%
return

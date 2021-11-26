
function [pop_data,pop_info] = compare_recordings(varargin)
% [pop_data,pop_info] = compare_recordings({filenames}, ... )
% Compare multiple ECAP or ENG recordings.
% 
% plots.compare_recordings uses persistent variables to minimise data
% loading overhead while exploring how best to visualise your data. 
% 
% Produces 1 page per active axon population and 1 row per configuration,
% roughly speaking. Much more useful for ENG than ECAP recordings.
% 
% filenames can be a cell array or a struct array as output by dir(). 
% 
% Input options:
% -data [pop_data, pop_info] : use pre-loaded data
% -clear          : clear persiste nt variables and return
% -opts {arglist} : argument list to underlying models.compose_waves
% -get-data       : return data and info without plotting
% 
% -pdf            : generate PDF comparison
% -one-row        : override info.group from models.compose_waves and show
%                   all of the spikerates/patterns on one row. 
% 
% V0.2 CDE - updated and added documentation and some options

named = @(v) strncmpi(v,varargin,length(v));
get_  = @(v) varargin{find(named(v))+1};

%% Parse input arguments - list of files to plot

if nargin == 0, list = {}; 
elseif isstruct(varargin{1}), list = varargin{1}; 
elseif iscell(varargin{1}), list = varargin{1}; 
else list = {};   
end

if isempty(list) && ~any(named('-data'))

  sel = tools.file('waves~/'); 
  list = {}; 
  while(ischar(sel))
    sel = uigetdir([sel filesep '..']);
    if ischar(sel), list(end+1) = {sel}; end %#ok<AGROW>
  end
  
  clear sel 
end

f_ = @(x) [x.folder filesep x.name];
if isstruct(list), list = arrayfun(f_,list,'unif',0); end

persistent data info

if any(named('-clear')), data = []; info = []; return
elseif any(named('-data')), data = get_('-data'); 
         info = varargin{find(named('-data'))+2};
         
elseif isempty(data) % Load data once
 
  if isempty(which('tools.make_SPARC_structure')), path(path,expdir('PN')), end
  
  data = []; 
  
  mcw_opts = {'-base','self'}; % options for models.compose_waves
  if any(named('-com')), mcw_opts = get_('-com'); 
  elseif any(named('-opt')), mcw_opts = get_('-opt'); 
  end
  if ~iscell(mcw_opts), mcw_opts = {mcw_opts}; end
  
  for ii = 1:numel(list)
    disp(['Reading ' strrep(list{ii},tools.file,'~')])
    [these,info] = models.compose_waves(list{ii}, mcw_opts{:});
    if isempty(data), data = these;
    else data = cat(3,data,these); 
    end
  end
  
  info.filename_list = list; 
end

if any(named('-get-data'))
  if nargout == 0, 
    assignin('caller','pop_data',data), 
    assignin('caller','pop_info',info),     
    clear, 
  end
  return
end

pop_data = data;
pop_info = info; 

clear sel ii f_ mcw_opts these 

%% Set up PDF and other options 

do.pdf = any(named('-pdf'));
do.plot_spectra = isfield(data,'hz');
do.plot_xc = ~do.plot_spectra;
do.zoom = []; 

if any(named('-zoom')), do.zoom = get_('-zoom'); end
if any(named('-one-row')), info.group(:) = 1; end
plots.PDF_tools('setup', do.pdf); 

pop_index = info.group;

%%

nA = size(data,1);
nRows = max(info.group); 
nY = size(data,3);
  
do.tweak = nRows < 8 && ~any(named('-no-tw')); 


colors = lines(max(7,size(data,3))); 
colors = [.5 .5 .5; colors];

for active_id = 1:nA  
  %% Load population-modulated spiking files

  axon_type = info.axon_type;
  active_type = info.active_ty{active_id} ; 
  
  time = sort(info.time);   
  get_roi_ = @(t) t >= time(1) & t <= time(end); 
  c_index = info.meta_col == 'c';
  if ~any(c_index), c_index(1) = true; end
  
  if do.plot_xc, title_str = 'Cross-correlation';
  else           title_str = 'Spectral';
  end
  
  %% Generate PDF pages
  clf
    
  for pp = 1:nRows
    %% Make each row on the image
  
    img = []; 
    coh = []; 
    iid = []; 
      
    for id = 1:nY % Build image and get average
      %%
      sel = (pop_index == pp);
      these = data(active_id,sel,id);       
      
      [c,order] = sort(info.metadata(sel,c_index));
      coh = [coh; c];               %#ok<AGROW>
      iid = [iid; 0*c+id];          %#ok<AGROW>

      c = colors(id,:);
        
      time = these(1).spk_time; % Go about constructing histogram
      roi = get_roi_(time); % units of ms

      px = time(roi); 
      py = mean(these(1).spk_rate(active_type,roi),1);
      % [px,py] = stairs_(px, py);
        dx = mean(diff(px))/2; 
        px = [1;1] * reshape(px,1,[]) + [-dx; dx];
        px = px([1 1:end end]);
        py = [1;1] * reshape(py,1,[]); 
        py = [0 py(:)' 0];
      % end stairs_

      subplot(nRows,3,3*pp - 2), hold on % plot histogram
      fill(px,py,c,'LineWidth',1.2,'EdgeColor',c,'FaceAlpha',0.4)

      subplot(nRows,3,3*pp), hold on % Plot average XC or spectrum

      if do.plot_xc
        px = these(1).lag;
        py = mean([these(order).xc],2);
        img = [img  these(order).xc]; %#ok<AGROW>
      elseif do.plot_spectra
        px = these(1).hz;
        py = mean([these(order).spect],2);
        img = [img  these(order).spect]; %#ok<AGROW>
      end
      plot(px,py,'LineWidth',1.3,'Color',c)
    end

    % Add annotations to spike-rate plot (left column)
    subplot(nRows,3,3*pp-2), axis tight, tools.tidyPlot
    if pp < nRows, set(gca,'XTickLabel',''); else xlabel('time (ms)'), end
    if pp == 1, title('Population spike-rate'), end
    if do.tweak, set(gca,'Position',get(gca,'Position') - [8 1 -8 -2]/100), end
    ylabel('imp/s/cell')

    % Add annotations to image plot (centre column)
    subplot(nRows,3,3*pp-1), cla reset      

    if do.plot_xc
      imagesc(these(1).lag,1:size(img,2),img')
      caxis([-1 1]*quantile(img(:),0.995))
    elseif do.plot_spectra        
      imagesc(these(1).hz,1:size(img,2),log10(img'))
      caxis(quantile(log10(img(:)),[0.005 0.995]))
    end
    hold on, axis tight xy, tools.tidyPlot
    if do.tweak, set(gca,'Position',get(gca,'Position') - [4 1 -7 -2]/100), end

    %% Clean up plots this row 
    
    if do.plot_xc % Set X-lim
      if ~isempty(do.zoom), xlim([-1 1]*do.zoom),
      else xlim(these(1).lag([1 end])) 
      end
      xlbl = 'lag (ms)';
    elseif do.plot_spectra, xlim([0 max(xlim)])
      xlbl = 'freq (Hz)';
    end
    if pp < nRows, set(gca,'XTickLabel',''); else xlabel(xlbl), end
    set(gca,'YColor','none'), xl = [min(xlim) range(xlim)/50];
    if pp == 1, title([title_str ' image']), end

    for id = 1:nY % add coherence indicator on y-axis
      c = colors(id,:);
      py = find(iid == id);
      px = xl(1) - xl(2)*log10(coh(iid==id)) - 1.5*xl(2);  
      axis(axis)
      % if do.plot_spectra, py = nB*py; end        
      plot(px(:), py(:), 'color',c,'LineWidth',1.3,'Clipping','off')
    end

    subplot(nRows,3,3*pp), axis tight, tools.tidyPlot
    if do.tweak, set(gca,'Position',get(gca,'Position') - [1 1 -7 -2]/100), end
    if pp < nRows, set(gca,'XTickLabel',''); else xlabel(xlbl),  end

    if do.plot_xc
      set(gca,'YAxisLocation','right')
      if ~isempty(do.zoom), xlim([-1 1]*do.zoom),
      else xlim(these(1).lag([1 end])) 
      end
    elseif do.plot_spectra
      set(gca,'YScale','log')
      xlim([min(xlim) xl(2)*50])
      if do.tweak, set(gca,'Position',get(gca,'Position') + [1 0 0 0]/100), end
    end
    if pp == 1, title([title_str ' average']), 
    elseif pp == nRows
    
      for id = 1:nY % add coherence indicator on y-axis
        c = colors(id,:);
        ylbl = regexp(info.filename_list{id},'(?<=\()[^\)]+','match','once');
        py = ylim*[1-id/nY; id/nY];         
        text(min(xlim),py,strrep(ylbl,'_','\_'),'color',c,'fontSize',8)        
      end

    end

  end % over population code freqs

  set(0,'ShowHiddenHandles','off')

  h = get(gcf,'Children');
  set(h(2:3:end),'YLim',[min([h(2:3:end).YLim]) max([h(2:3:end).YLim])])

  txt = sprintf(', %s', axon_type{active_type});  
  txt = sprintf('%s: Population spike-rate variation', txt(3:end));
  tools.suptitle(txt)

  %%
  if do.pdf, plots.PDF_tools(gcf, 'a%03d-summary%d.ps', active_id, 0); 
       pause(0.05), close(gcf)
  else pause(0.05), figure
  end
end % active_ID

if do.pdf
  plots.PDF_tools('combine', 'Compare_population_responses (%d).pdf');  
else close(gcf), 
  if ~isempty(which('tileMyFigures')), tileMyFigures, end
end

if nargout == 0
  assignin('caller','results',pop_data)
  assignin('caller','pop_info',pop_info)
  clear, return
end


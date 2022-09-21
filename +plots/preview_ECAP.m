

function preview_ECAP(varargin)
% ECAP_response(waves_folder, ...)
% 
% -pdf : make a pdf
% -eidors-file : use specified EIDORS file for context 
%               (default: read INPUTS file)
% -raster : make ECAP raster (old default view)
% 
% The default plots.response_wave isn't really what we want to look at ECAP
% responses, so this function provides two better views. The first set of
% figures was generated in the context of the 4x2 array experiments, and
% looks at the response on a given pair across spike-rates for a couple
% bipolar pairs per page. 
% 
% The second, accessed using -raster, looks at every 2-electrode pair but
% shows only a single spike-rate per page. the second view extracts spectra
% and tries to look at cross-fasciclar sensitivity isn't stellar. also
% shows underling raster plot of the simulated ECAP response.
% 
% v0.2 25-Jun-2020 CDE
% 

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('-raster')) % OLD default view for this  
  make_ECAP_raster(varargin{:})  
  return
end

do_PDF = any(named('-pdf'));
plots.PDF_tools('setup',do_PDF);

%%

[ecap_folder,EM,nerve] = tools.parse_arguments(varargin,'wave+context','waves~/','LOAD');


nE = sum(strncmpi({EM.model.electrode.name},'Elec',4)); 
nF = sum(strncmpi(fieldnames(EM),'Fascicle',4));

elec_image = [1 2; 3 4; 1 3; 2 4]; 
if any(named('-diag')), elec_image = [elec_image; 2 3; 1 4]; end
if any(named('-elec')), elec_image = get_('-elec'); end

elec_map = elec_image; 

while max(elec_map) < nE  
  elec_map = [elec_map; elec_image + max(elec_map(:))]; %#ok<AGROW>
end

ok = all(elec_map <= nE,2);
elec_map = elec_map(ok,:); 

%% Load in all permutations of the data (kinda tedious, avoid if possible)

if evalin('caller','exist(''ALL_response'',''var'')')
  
  ALL_response = evalin('caller','ALL_response');
  FASC_response = evalin('caller','FASC_response');
  pop_info = evalin('caller','pop_info');
else

  % Calling plots.response_waves( ..., -gather-data ) uses that machinery
  % to load all of the relevent data for us and preprocess it correctly. 
  % 
  % This call generates 2 variables, pop_info / pop_response, which would
  % usually get used to plot the response waves. 
  
  mcw_args = {'-chan',elec_map,'-quick'};
  if any(named('-base')), mcw_args = [mcw_args {'-base', get_('-base')}]; end
  if any(named('-args')), more_args = get_('-args'); % generic passthrough
    if iscell(more_args), mcw_args = [mcw_args more_args];
    else                  mcw_args = [mcw_args {more_args}];
    end
  end
  
  FASC_response = cell(nF,1);
  
  for ff = 1:nF
   [FASC_response{ff},~] = models.compose_waves(ecap_folder,mcw_args{:},'-fasc',ff);
  end  
 [ALL_response, pop_info] = models.compose_waves(ecap_folder,mcw_args{:});
  pop_info.color = cat(1,EM.axons.pop.color);
  
  assignin('caller','ALL_response',ALL_response)
  assignin('caller','FASC_response',FASC_response)
  assignin('caller','pop_info',pop_info)
  clc, disp('all data loaded')
end

%%

% the output PDF has 2 sets of pages -
% Group A) for each 4 electrodes, make a 1-page summary showing how the 
% responses change as a function of spike-rate for that set of electrodes. 

nRows = size(elec_image,1);
if any(named('-rows')),  nRows = get_('-rows'); end
nPages = ceil(size(elec_map,1) / nRows);

opts.recording_map = elec_map; 
opts.fascicle = 1:nF;
opts.t_roi = [-15 50];
opts.a_roi = [1.5 20];
if any(named('-time')), opts.t_roi = get_('-time'); end
if any(named('-roi')),  opts.a_roi = get_('-roi');  end

opts.chan = 1;

for page_id = 1:nPages
  
  figure, clf
  set(gcf,'Position',[745 90 770 960])
  
  for row_id = 1:nRows
    
    opts.chan = (page_id-1)*nRows + row_id;   
    if opts.chan > size(elec_map,1), break, end
    subplot(nRows,1,row_id)
    make_response_row(EM,nerve,ALL_response,pop_info,opts);
  end
  
  %% apply formatting
  
  h = findobj(gcf,'UserData','wave');
  p = cat(1,h.Position);
  set(h,'YLim',[min([h.YLim]) max([h.YLim])])
  for ii = 1:length(h)    
    if p(ii,2) == min(p(:,2)), continue, end
    set(h(ii),'XTickLabel','','YTickLAbel','')
    xlabel(h(ii),''), ylabel(h(ii),'')
  end
  linkaxes(h,'x')
  
  h = findobj(gcf,'UserData','errorbar');
  p = cat(1,h.Position); y = cat(1,h.YLim);
  set(h,'YLim',[min([h.YLim]) max([h.YLim])])
  for ii = 1:length(h)    
    
    if y(ii,2) ~= max(y(:,2))      
      delete(findobj(h(ii),'type','text'))
    end    
    if p(ii,2) == min(p(:,2)), continue, end
    set(h(ii),'XTickLabel','','YTickLAbel','')
    xlabel(h(ii),''), ylabel(h(ii),'')
  end
  
  h = findobj(gcf,'UserData','nerve');
  p = cat(1,h.Position);
  set(h,'XtickLabel','');
  
  for ii = 1:length(h)    
    if p(ii,2) == max(p(:,2)), ylabel(h(ii),'(mm)'), continue, end
    set(h(ii),'YTickLAbel','')    
  end
  
  idx = (page_id-1)*nRows + (1:nRows); 
  idx(idx > size(elec_map,1)) = [];
  erange = elec_map(idx,:);
  suptitle(sprintf('E%d-E%d Bipolar Recordings', min(erange(:)),max(erange(:))));
  
  %% make PDF page
  plots.PDF_tools(gcf,do_PDF,'A-%03d-waves.ps',page_id,'-move','-resize','-portrait'); 
  
end

%%

% Group B) for each electrode config, make a page showing how the response
% changes for each fascicle in the recording set. 

for page_id = 1:size(elec_map,1)
  
  figure, clf
  set(gcf,'Position',[745 90 770 960])
  
  for ff = 1:nF
    
    opts.chan = page_id;       
    opts.fascicle = ff; 
    
    subplot(nF,1,ff)
    make_response_row(EM,nerve,FASC_response{ff},pop_info,opts);
  end
  
  %% apply formatting
  
  h = findobj(gcf,'UserData','wave');
  p = cat(1,h.Position);
  set(h,'YLim',[min([h.YLim]) max([h.YLim])])
  for ii = 1:length(h)    
    if p(ii,2) == min(p(:,2)), continue, end
    set(h(ii),'XTickLabel','','YTickLAbel','')
    xlabel(h(ii),''), ylabel(h(ii),'')
  end
  linkaxes(h,'x')

  h = findobj(gcf,'UserData','errorbar');
  p = cat(1,h.Position); y = cat(1,h.YLim);
  % set(h,'YLim',[min([h.YLim]) max([h.YLim])])
  for ii = 1:length(h)    
    
    if y(ii,2) ~= max(y(:,2))      
      delete(findobj(h(ii),'type','text'))
    end    
    if p(ii,2) == min(p(:,2)), continue, end
    set(h(ii),'XTickLabel',''); % ,'YTickLAbel','')
    xlabel(h(ii),''), ylabel(h(ii),'')
  end
  
  h = findobj(gcf,'UserData','nerve');
  p = cat(1,h.Position);
  set(h,'XtickLabel','');
  
  for ii = 1:length(h)    
    if p(ii,2) == max(p(:,2)), ylabel(h(ii),'(mm)'), continue, end
    set(h(ii),'YTickLAbel','')    
  end
  
  suptitle(sprintf('E%d.%d Fascicle sensitivity', elec_map(page_id,:)))
  
  %%
  plots.PDF_tools(gcf,do_PDF,'B-%03d-waves.ps',page_id,'-move','-resize','-portrait'); 
end

plots.PDF_tools('compile',do_PDF,'RMS Response Curves (%d).pdf')

%%
return



function make_ECAP_raster(varargin)

% note : ECAP_response -raster is completely different, and uses the
% following options: 

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};


do_PDF = any(named('-pdf'));
f_ = @(x) [x.folder filesep x.name];

plots.PDF_tools('setup',do_PDF);


% Get primary data

if any(named('-dir')), data_ = get_('-dir'); % user specified folder
  if any(data_ == '~'), data_ = tools.file(data_); end
  if ~ismember(data_(end),filesep), data_ = [data_ filesep]; end
elseif nargin > 0 && varargin{1}(1) ~= '-' % leading entry
  data_ = varargin{1}; 
  wp = tools.file('sub~\waves\');
  if any(data_ == '~'), data_ = tools.file(data_);
  elseif ~isfolder(data_) && isfolder([wp data_]), data_ = [wp data_]; 
  end  
  if ~ismember(data_(end),filesep), data_ = [data_ filesep]; end
else
  data_ = dir(tools.file('sub~\waves\ecap*'));
  if numel(data_) == 1, data_ = [f_(data_) filesep];
  else data_ = '';
  end  
end

if ~isfolder(data_)  
  [data_] = uigetdir(tools.file('sub~\waves\'));
  if ~ischar(data_), return, end % Selecion cancelled, exit gracefully
  if ~ismember(data_(end),filesep), data_ = [data_ filesep]; end
end

assert(contains(data_,'stim'),'Please select a .../waves/stim folder'); 

%% Get Support files

list = dir(tools.file('sub~\eidors\sensitivity*.mat'));
if any(named('-eidors-file')), eidors_file = get_('-eidors-file');
else eidors_file = tools.INPUT_file(list, data_);  
end

EM = load(eidors_file); 

for f = fieldnames(EM.utils)' % Create utility local functions
    EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
    eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
end

f = dir(tools.file('sub~/axons/axons*.mat')); 
axons_file = tools.INPUT_file(f, eidors_file); 
EM.axons = load(axons_file);

f = dir(tools.file('~/source/fascicles/*.splines.dat')); 
spline_file = tools.INPUT_file(f, axons_file); 

nerve = read_dat_file(spline_file,'index',EM.info.SplineIndex);
nerve.fascicle = nerve.outline / 1000; % convert to mm
nerve.fascicle(:,2,:) = nerve.fascicle(:,2,:) + 0.050; % add offset

data_ = @(f) [data_ f];

%% Get File list 

list = dir(data_('*.mat'));
rates = regexp({list.name},'_k[\d\.]+_','match','once'); 
list(cellfun(@isempty,rates)) = []; 
rates(cellfun(@isempty,rates)) = []; 

[rates,~,rate_id] = unique(rates);
[rates,order] = sort(cellfun(@(s) str2double(s(3:end-1)), rates));
order(order) = 1:numel(order);
rate_id = order(rate_id);

[rate_id,order] = sort(rate_id);
rate_id = rate_id';
list = list(order);
clear order

%%

chronux_opts = tools.setupChronux; 
chronux_opts.fpass = [0.2 3e3];

axon_color = EM.axons.pop.axon_color; %#ok<NASGU>
n_elec = numel(EM.model.electrode)-1;

style = {'linewidth',1.1,'color'};
C = lines(7); G = @(v) [v v v]/10; %#ok<NASGU>

t_roi = [-15 50];

for r_id = 1:length(rates)
  
  this_rate = find(rate_id == r_id)';
  D = load(data_(list(this_rate(1)).name));
  
  %% Raster plot (untyped)
  
  figure, clf
  subplot(3,3,[2 3]), cla, hold on

  z0 = 0; 
  for ii = 1:numel(D.axontype)
    ok = (D.raster{ii}.spk_time > t_roi(1) & D.raster{ii}.spk_time < t_roi(end));
    plot(D.raster{ii}.spk_time(ok), D.raster{ii}.spk_axon(ok) + z0,'.', ...
                'Color',G(2))
              % 'Color',axon_color(ii,:))
    z0 = z0 + length(D.raster{ii}.axon_group);
  end

  axis tight, xlim(t_roi), tools.tidyPlot, 
  set(gca,'YTick',[],'TickLength',[1 1]/150)
  ylabel(sprintf('%0.1f imp/s', rates(r_id)))
  
  %% Response waves

  chan_color = jet(n_elec) .* [1 .85 1];  
  % dy = quantile(D.waves(:),0.8); % µv
  
  subplot(3,3,[5 6]), cla, hold on
  ok = (D.time > t_roi(1)& D.time < t_roi(end));
  
  for ee = 1:2:n_elec
    plot(D.time(ok),sum(D.waves(ok,ee,:) - D.waves(ok,ee+1,:),3), ...
              style{:},chan_color(ee,:))
  end
  
  % plot(D.time(ok),sum(D.waves(ok,3,:) - D.waves(ok,4,:),3),style{:},G(5))
  % plot(D.time(ok),sum(D.waves(ok,1,:) - D.waves(ok,2,:),3)-dy,style{:},G(2))
  % for ii = 1:numel(D.axontype)
  %     plot(D.time(ok),D.waves(ok,1,ii)-D.waves(ok,2,ii)+dy*ii, ...
  %                               style{:},axon_color(ii,:))
  %   end
  
  xlim(t_roi), 
  tools.tidyPlot, set(gca,'TickLength',[1 1]/150)
  xlabel('time (ms)'), ylabel('(µV)')
  
  % SCALEBAR
  % set(gca,'YColor','none')
  % plot([1 1]*min(xlim),[0 -2],style{:},get(gca,'XColor'))
  % text(xlim*[1.01;-0.01],-1,'2 µV','Color',get(gca,'XColor'),'Rotation',90, ...
  %      'horizontalalignment','center','verticalalignment','bottom','clipping','off')

  %% Array Drawing  
  subplot(4,3,1), cla reset, hold on
  sel = EM.model.object_id{strncmpi(EM.model.object_name,'PDMS',4)}; 
  idx = EM.model.elems(sel,:);
  idx = idx(sum(reshape(y_(idx),[],4) >= -10*eps,2) > 2,:);

  patch('Faces',idx, 'Vertices',EM.model.nodes(:,[1 3]), ...
        'EdgeColor',[.7 .7 .7], 'FaceColor','w','EdgeAlpha',0.2)

  outline = EM.model.nodes(unique(idx),[1 3]);
  idx = convhull(outline); 

  plot(outline(idx,1),outline(idx,2), '-', style{:},G(5))

  for ee = 1:n_elec
    
    e_n = EM.model.electrode(ee).nodes;
    e_xy = [min(x_(e_n)) min(z_(e_n)) range(x_(e_n)) range(z_(e_n))]; 
    
    if mod(ee,2) == 0, fc = (chan_color(ee,:) + 0.2)/1.2;
    else               fc = (chan_color(ee,:) + 1.2)/2.2;
    end
    rectangle('Position',e_xy,'FaceColor',fc,'EdgeColor', ...
                chan_color(ee,:),'LineWidth',1.2)
  end
  
  axis equal, axis off, axis(axis * [1.05 -0.05 0 0; -0.05 1.05 0 0; ...
                                     0 0 1.05 -0.05; 0 0 -0.05 1.05]);
  for ff = 1:size(nerve.fascicle,3)    
    
    
    e_xy = [min(xlim) min(nerve.fascicle(:,1,ff)) ...
            max(xlim) max(nerve.fascicle(:,1,ff))]; 
    
    fill(e_xy([1 3 3 1 1]),e_xy([2 2 4 4 2]),G(2),'EdgeColor', ...
                G(2),'LineWidth',0.8,'FaceAlpha',0.3)
  end
  % ylim(ylim - 0.25*range(ylim))
  
  %% Spectra (by electrode)

  chronux_opts.Fs = 1000 / mean(diff(D.time)); % units of kS/s

  avg_spectrum = []; 
  % alpha = 0.1;

  subplot(3,3,[8 9]), cla, hold on

  for rr = this_rate    
    D = load(data_(list(rr).name));      
    assert(size(D.waves,2) == n_elec,'%s: Electrode count mismatch', list(rr).name)
    
    for ee = 1:2:n_elec, idx = (ee+1)/2;   
        
      % wave = [squeeze(D.waves(:,ee-1,:)-D.waves(:,ee,:))     ...  
      %            sum(D.waves(:,ee-1,:)-D.waves(:,ee,:),3)];
        
      wave = squeeze(sum(D.waves(:,ee,:,:)-D.waves(:,ee+1,:,:),4));
      wave = [wave sum(wave,2)]; %#ok<AGROW>
      [spec,hz] = mtspectrumc(wave, chronux_opts);

      if isempty(avg_spectrum), avg_spectrum = spec;
           avg_spectrum(end,1,n_elec/2) = 0;
      else avg_spectrum(:,:,idx) = avg_spectrum(:,:,idx) + spec; %#ok<AGROW>
      end
    end
  end

  avg_spectrum = avg_spectrum / rr; 

  for ee = 1:n_elec/2,
    plot(hz,avg_spectrum(:,end,ee),style{:},chan_color(2*ee-1,:))    
  end
  
  % for ii = 1:numel(D.axontype)
  %   plot(hz,avg_spectrum(:,ii),style{:},axon_color(ii,:))
  % end
  % plot(hz,avg_spectrum(:,end),style{:},G(2))

  xlabel('frequency, Hz'), ylabel('Spectral power')
  set(gca,'YScale','log'), tools.tidyPlot
  ylim(10.^[-9.5 0])
  
  %% Fascicle sensitivity plots 
  
  n_fasic = size(nerve.fascicle,3);
  
  subplot(4,3,4), cla, hold on
  for ff = 1:n_fasic
    fill(nerve.fascicle(:,1,ff), nerve.fascicle(:,2,ff), G(8), ... 
          'LineWidth',1.2,'EdgeColor',G(2))
    text(mean(nerve.fascicle(:,1,ff)), mean(nerve.fascicle(:,2,ff)), ... 
          sprintf('f#%d',ff),'Color','k','Horiz','center','fontsize',8)
  end
  axis image, tools.tidyPlot, grid on
  if min(ylim) > 0, ylim([-0.03 1] * max(ylim)), end
  set(gca,'XTick',-1:0.1:1,'YTick',-1:0.1:1)

  subplot(4,3,[7 10]), cla, hold on
  
  
  roi_A = (hz > 600);
  roi_C = (hz > 10 & hz < 300);
  
  est_C = squeeze(sum(avg_spectrum(roi_C,:,:),1));
  est_A = squeeze(sum(avg_spectrum(roi_A,:,:),1));
  
  for ff = 1:n_fasic
    idx = convhull(est_C(ff,:),est_A(ff,:));
    fill(est_C(ff,idx),est_A(ff,idx),G(5),'EdgeColor',G(6),'FaceAlpha',0.2); 
  end
  
  for ee = 1:n_elec/2,
    plot(est_C(1:end-1,ee),est_A(1:end-1,ee),'o',style{:},chan_color(2*ee-1,:))    
  end
  
  for ff = 1:n_fasic    
    text(max(est_C(ff,:)), mean(est_A(ff,:)), ... 
          sprintf('f#%d',ff),'Color','k','fontsize',12)
  end
  
  set(gca,'xscale','log','yscale','log')
  tools.tidyPlot
  xlabel('C-fibre (10-300 hz)')
  ylabel('A-fibre (600+ hz)')
  
  plots.PDF_tools(gcf,do_PDF,'waves-%03d.ps',r_id);
  
end


if do_PDF
  plots.PDF_tools('compile',[eidors_file '-ecaps (%d).pdf'])
end






% Each row of the (default) ECAP_response PDF
function make_response_row(EM,nerve,results,info,opts)

p = get(gca,'Position') - [0.05 0 -0.05 0]; delete(gca);

nE = sum(strncmpi({EM.model.electrode.name},'Elec',4)); 

electrode_color = jet(nE) .* [1 .85 1];  
C = lines(7); G = @(v) [v v v]/10;  %#ok<NASGU>

for f = fieldnames(EM.utils)' % Create utility local functions
    EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
    eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
end

%% Fascicle position plots 
axes('Position',p .* [1 1 0.25 0.45] + [0 p(4)*0.65 0 0]), cla, hold on


if range(nerve.outline(:)) > 2*max(EM.info.DomainSize)
  nerve.outline = nerve.outline / 1e3; % um to mm
end

nF = size(nerve.outline,3);
for ff = 1:nF
  if ismember(ff,opts.fascicle), 
       style = {G(8),  'LineWidth',1.2,'EdgeColor',G(2)};
  else style = {G(9.5),'LineWidth',1.1,'EdgeColor',G(6)};
  end  
  fill(nerve.outline(:,1,ff), nerve.outline(:,2,ff), style{:})
end

axis image, tools.tidyPlot, grid on
if min(ylim) > 0, ylim([-0.03 1] * max(ylim)), end
set(gca,'XTick',-1:0.1:1,'YTick',-1:0.1:1)
set(gca,'UserData','nerve')  


%% Electrode image
axes('Position',p .* [1 1 0.25 0.5]), cla, hold on
style = {'LineWidth',1.2,'Color'};

sel = EM.model.object_id{strncmpi(EM.model.object_name,'PDMS',4)}; 
idx = EM.model.elems(sel,:);
idx = idx(sum(reshape(y_(idx),[],4) >= -10*eps,2) > 2,:);

patch('Faces',idx, 'Vertices',EM.model.nodes(:,[1 3]), ...
      'EdgeColor',[.7 .7 .7], 'FaceColor','w','EdgeAlpha',0.2)

outline = EM.model.nodes(unique(idx),[1 3]);
idx = convhull(outline); 
plot(outline(idx,1),outline(idx,2), '-', style{:},G(5))

for ee = 1:nE

  e_n = EM.model.electrode(ee).nodes;
  e_xy = [min(x_(e_n)) min(z_(e_n)) range(x_(e_n)) range(z_(e_n))]; 

  fc = electrode_color(ee,:); 
  ec = electrode_color(ee,:); 

  if opts.recording_map(opts.chan,1) == ee,  fc = (fc + 0.3) / 1.3;
  elseif opts.recording_map(opts.chan,2) == ee, fc = (fc + 1.5) / 2.5; 
  else fc = (fc + 3)/5; ec = (ec + 2) / 4;
  end

  rectangle('Position',e_xy,'FaceColor',fc,'EdgeColor', ec, ...
                            'LineWidth',1.2)
end

axis equal, axis off, axis(axis * [1.05 -0.05 0 0; -0.05 1.05 0 0; ...
                                   0 0 1.05 -0.05; 0 0 -0.05 1.05]);
for ff = 1:size(nerve.outline,3)    
  e_xy = [min(xlim) min(nerve.outline(:,1,ff)) ...
          max(xlim) max(nerve.outline(:,1,ff))]; 

  fill(e_xy([1 3 3 1 1]),e_xy([2 2 4 4 2]),G(2),'EdgeColor', ...
              G(2),'LineWidth',0.8,'FaceAlpha',0.3)
end
set(gca,'UserData','array')
clear ee ff e_xy ec fc idx sel 

%%

axes('Position',p .* [1 1 0.38 1.05] + [0.26 0 0 0]), cla, hold on

color = tools.magma(ceil(size(info.spikerate,1)*1.1) + 2); 

fs = 1./mean(diff(info.time));
time = (1:size(results(1).wave,1)) / fs; 
time = time - mean(time(:));
roi = (time > min(opts.t_roi) & time <= max(opts.t_roi));

wave = arrayfun(@(x) x.wave(roi,opts.chan), results, 'unif',0);

for gg = 1:max(info.group)
  
  wsel = [wave{end,info.group == gg}];
  
  y0 = mean(wsel,2)'; 
  ye = std(wsel,[],2)'; 
  
  plot(time(roi),y0,'Color',color(gg,:),'LineWidth',1.2);
  fill([time(roi) fliplr(time(roi))], [y0+ye fliplr(y0-ye)], color(gg,:), ...
      'FaceAlpha',0.2,'EdgeColor','none');
end

axis tight, tools.tidyPlot
set(gca,'UserData','wave')

rectangle('Position',[min(opts.a_roi) min(ylim) range(opts.a_roi) range(ylim)], ...
          'FaceColor',G(9),'edgeColor','none')
  
set(gca,'Children',flipud(get(gca,'Children')),'layer','top')
ylabel('µV'), xlabel('ms')

%%

roi = (time > min(opts.a_roi) & time <= max(opts.a_roi));
wave = arrayfun(@(x) x.wave(roi,opts.chan), results, 'unif',0);
y0 = zeros(max(info.group),size(results,1));
ye = zeros(max(info.group),size(results,1));

if isfield(info,'stimulus')
  x0 = info.stimulus.current; xlbl = 'Stimulus (µA)';
  
else x0 = info.spikerate; xlbl = 'Spike-rate';
end

for gg = 1:max(info.group)
  
  for tt = 1:size(results,1)
  
    wsel = [wave{tt,info.group == gg}];
    wrms = sqrt(mean(wsel.^2,1));

    y0(gg,tt) = mean(wrms); 
    ye(gg,tt) = std(wrms);
  end
end

axes('Position',p .* [1 1 0.27 1.05] + [0.63 0 0 0]), cla, hold on

for tt = 1:size(results,1)
  
  if tt > size(info.color,1), color = G(2); lbl = 'All';
  else color = info.color(tt,:); lbl = info.axon_type{tt};
  end
  
  errorbar(x0,y0(:,tt),ye(:,tt),'Color',color,...
      'LineWidth',1.2,'CapSize',4,'marker','o','markerFAceColor',color,'markersize',5)
  text(info.spikerate(end)+3,y0(end,tt),lbl,'Color',color,'FontSize',7)
end

tools.tidyPlot, xlabel(xlbl), ylabel('V_{RMS} (µV)')
set(gca,'UserData','errorbar')


%%
return

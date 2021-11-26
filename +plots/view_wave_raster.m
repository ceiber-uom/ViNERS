
function view_wave_raster(varargin)
% function view_wave_raster(filename, varargin)
% 
% This function generates the invididual panels for EMBC fig. 3
% 
% Options: 
%  filename = '?' : use UI file picker (default if nargin = 0)
% -ax [axon_file] : use specified axon file (default: newest)
% -roi [xlim]     : set ROI for display
% -fascicle [#]   : filter fascicles (default: sum over all)
% -dy [delta µV]  : set Y offset for display
% 
% V0.2 CDE 18 Nov 2021


named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

%% Figure illustrating responses [Single spike-rate example]

if nargin > 0, data_file = varargin{1}; else data_file = ''; end
if contains(data_file,'~'), data_file = tools.file(data_file); end
if ~exist(data_file,'file')
    data_file = tools.file('sub~/waves/*.mat','-prompt');
    if ~exist(data_file,'file'), return, end % cancelled
end

[data_dir,data_file] = fileparts(data_file); 
data_dir = [data_dir filesep];

D = load([data_dir data_file '.mat']);

ax_file = tools.file('get','sub~/axons/axon*.mat','newest');
if any(named('-ax')), ax_file = get_('-ax'); end
load(ax_file,'pop');

D.axontype = D.options.class;
D.axon_color = cat(1,pop.color);
% TODO - set electrode grid, plot opts for raster, types
if any(named('-pop-col')), D.axon_color = get_('-pop-col'); end


if iscell(D.raster), D.raster = [D.raster{:}]; end

%%

t_roi = max(abs(D.time)) - 25; 
if any(named('-roi')), t_roi = get_('-roi');
  if ischar(t_roi), t_roi = str2double(t_roi); end
  if all(isnan(t_roi)), t_roi = 25; end
end
if numel(t_roi) == 1, t_roi = [-1 1]*abs(t_roi); end

clf, subplot(3,1,1), cla, hold on

z0 = 0; 
for ii = 1:numel(D.axontype)
  
    if isstruct(D.raster) && isfield(D.raster,'spike') % from models.ecap_recording
      spk_time = {D.raster(ii).spike.init}; 
      spk_axon = arrayfun(@(n) n*ones(size(spk_time{n})), 1:numel(spk_time),'unif',0);      
      spk_time = cat(1,spk_time{:}); spk_axon = cat(1,spk_axon{:}); 
      if isempty(spk_time), continue, end,       
      spk_time(:,2) = []; spk_axon(:,2) = []; 
    else
      spk_time = D.raster(ii).spk_time;
      spk_axon = D.raster(ii).spk_axon;
    end
    
    ok = (spk_time >= t_roi(1) & spk_time < t_roi(end));

    if any(named('-raster-o')), 
      axon_y = get('-raster-o'); 
      z_increment = nanmax(axon_y);
      axon_y = axon_y(spk_axon(ok));       
    else axon_y = spk_axon(ok);
      z_increment = length(D.raster(ii).axon_group);
    end
    
    plot(spk_time(ok), axon_y + z0,'.','Color', D.axon_color(ii,:))
    z0 = z0 + z_increment;
end

tools.tidyPlot, set(gca,'YTick',[],'TickLength',[1 1]/150), xlim(t_roi)

% Row 2 - response waves

dy = quantile(abs(D.waves(abs(D.waves)  > 1e-12)),0.99); % µv

if any(named('-dy')), dy = get_('-dy'); end

if size(D.waves,4) > 1
   
  if any(named('-f')), fid = get_('-f'); else fid = 1:size(D.waves,3); end
  D.waves = permute(sum(D.waves(:,:,fid,:),3),[1 2 4 3]);
end


elec = [1 2; 3 4]; 

if any(named('-elec')), elec = get_('-elec');
elseif any(named('-chan')), elec = get_('-chan');
elseif any(named('-pair')), elec = get_('-pair');
end
 
  

subplot(3,1,2), cla, hold on
style = {'linewidth',1.1,'color'};
C = lines(7); G = @(v) [v v v]/10; %#ok<NASGU>

ok = (D.time > t_roi(1)& D.time < t_roi(end));

if size(elec,2) > 1
     get_wave_ = @(w,e) sum(w(:,e(1),:) - mean( ...
                            w(:,e(e>0 & (1:numel(e))>1),:),2),3);
else get_wave_ = @(w,e) sum(w(:,e(1),:), 3); 
end
                     
if size(elec,1) > 1
  for ee = 2:size(elec,1)
    plot(D.time(ok),(ee-1)*dy + ...
         get_wave_(D.waves(ok,:,:), elec(ee,:)), style{:},G(5))
  end
else
  
  
end
plot(D.time(ok),get_wave_(D.waves(ok,:,:), elec(1,:)),style{:},G(2))

for pp = 1:numel(D.axontype)
  plot(D.time(ok),get_wave_(D.waves(ok,:,pp),elec(1,:)) - dy*pp, ...
       style{:},D.axon_color(pp,:))  
end

axis tight, tools.tidyPlot, xlim(t_roi)
set(gca,'YColor','none','TickLength',[1 1]/150)

sb = min(10.^round(log10(max(-ylim))), 10);

plot((xlim*[1.01;-0.01])*[1 1],[0 sb]+min(ylim),style{:},get(gca,'XColor'),'Clipping','off')
text(xlim*[1.025;-0.025],sb/2+min(ylim),sprintf('%g µV', sb),'Color',get(gca,'XColor'),'Rotation',90, ...
    'horizontalalignment','center','verticalalignment','bottom')

xlabel('time (ms)')  
linkaxes(get(gcf,'children'),'x')
  
%%
subplot(3,1,3), cla reset, hold on

chronux_opts = tools.setupChronux; 
chronux_opts.Fs = 1000 / mean(diff(D.time)); % units of kS/s
chronux_opts.fpass = [0.2 4e3];

avg_spec = []; 
data_list = dir([data_dir regexprep(data_file,'k\d+','*')]);

if isempty(data_list)  
  data_list = dir([data_dir regexprep(data_file,'\(\d+\)','*')]);
end

if numel(data_list) <= 1
  warning('ViNERS:viewWaves:cannotAverage',... 
        'Only %d files found which are replicates of "%s", cannot produce an average.', ... 
        numel(data_list), data_file)

  all_waves = get_wave_(D.waves(:,:,:),elec(1,:));   
  for pp = numel(D.axontype):-1:1
    all_waves = [get_wave_(D.waves(:,:,pp),elec(1,:)) all_waves]; %#ok<AGROW>
  end
      
  if isempty(which('mtspectrumc'))
   [avg_spec,hz] = pwelch(all_waves,[],[],[],chronux_opts.Fs);
    avg_spec(hz > hz(end)/2,:) = []; 
    hz(hz > hz(end)/2,:) = []; 
  else [avg_spec,hz] = mtspectrumc(all_waves,chronux_opts);
  end
end

%%
for ff = 1:length(data_list)  
  
  D = load([data_dir data_list(ff).name]);
 
  D.axontype = D.options.class;
  D.axon_color = cat(1,pop.color);
  

  if size(D.waves,4) > 1

    if any(named('-f')), fid = get_('-f'); else fid = 1:size(D.waves,3); end
    D.waves = permute(sum(D.waves(:,:,fid,:),3),[1 2 4 3]);

  end

  if strcmp(data_list(ff).name,data_file), alpha = 0.6;
  else                                     alpha = 0.1;
  end
  

  all_waves = get_wave_(D.waves,elec(1,:));   
  for pp = numel(D.axontype):-1:1
    all_waves = [get_wave_(D.waves(:,:,pp),elec(1,:)) all_waves];
  end  
  
  if isempty(which('mtspectrumc'))
      [spec,hz] = pwelch(all_waves,[],[],[],chronux_opts.Fs);
      spec(hz > hz(end)/2,:) = []; 
       hz(hz > hz(end)/2,:) = []; 
  else [spec,hz] = mtspectrumc(all_waves, chronux_opts);
  end

  for ii = 1:numel(D.axontype)
    plot(hz,spec(:,ii),'color',[D.axon_color(ii,:) alpha])
  end
  % plot(hz,spec(:,end),style{:},G(2))

  if isempty(avg_spec), avg_spec = spec;
  else avg_spec = avg_spec + spec;
  end
end

if ~isempty(data_list), avg_spec = avg_spec / ff; end

for ii = 1:numel(D.axontype)
  plot(hz,avg_spec(:,ii),style{:},D.axon_color(ii,:))
end
plot(hz,avg_spec(:,end),style{:},G(2))

xlabel('frequency, Hz'), ylabel('Spectral power')
set(gca,'YScale','log'), tools.tidyPlot

% ylim(10.^[-9.5 0])

tools.suptitle(strrep(data_file,'_','\_'))

  
%%

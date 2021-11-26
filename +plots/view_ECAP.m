
function data = view_ECAP(varargin)
% ECAP_data  = view_ECAP(filename, ...)
% 
% Options: 
% -data [data]       : use previously loaded data 
%                      (e.g. with different filtering)
% -opts {list}       : pass the specified list as options to load_ECAPs
%                     (see below)
% -index [#]         : load the #th ECAP folder (alphabetical)
% -noise [rms]       : amount of noise, default 0.275 µV RMS pre filtering
% -no-filter         : suppress default filtering 
% -sig               : compute sigmoid (logistic) curve
% -roi [4.5 9]       : set analysis ROI (ms post stimulus)
% -load              : invoke load_ECAPs and pass data to caller
% 
% DISPLAY options 
% -view [0 12]       : set display ROI (ms post stimulus)
% -y-step [1]        : set y step for ECAP display (default 1 µV)
% -xl [80 2000]      : set xlim for response curve (µA stimulus)
% -extend            : plot un-clipped response curve 
% 
% ==== load_ECAPs (internal function) ====
% DATA SELECTION options
% -dir [ecap_path]   : use specified ECAP path (default: "waves~")
% -?, ?              : use a UI menu based on the content of ecap_path
% -index [#]         : load the #th ECAP folder (alphabetical)
% -all               : load ALL ECAP data folders 
% -fast              : load summation of all axon classes (1:4) only 
%                      (default: load 1,2,3,4 and [1:4])
% -by-f              : load data for each fascicle seperately 
% 
% PRE-PROCESSING option
% -filter              : apply bandpass or other filter to wave
% -hz [low high order] : apply bandpass filtering (also supports lowpass,
%                        highpass, and variable-order)
%                        default [500 3000 2] Hz
% -zerolag             : apply zero-lag fitering (default: time-reversed,
%                        set negative order for causal filtering)
% -noise               : set noise amount 
%                        (default 0.0275 µV or 0.275 before filtering)
% -debug-filter        : visualisation of filtering 
% -args {arglist}      : pass additional arguments to models.compose_waves
% 
% MISC optins
% -q                 : quiet mode
% --d [data]         : reprocess specified data
% 
% last updated 9 Nov 2021 Calvin Eiber
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

%% Load the data 

if any(named('-d')), data = get_('-d'); 
elseif ~exist('data','var')
  
  if any(named('-opts')), load_args = get_('-opts'); 
    if ~iscell(load_args), load_args = {load_args}; end
  elseif any(named('-index')), load_args = {'-index', get_('-index')};  
  elseif nargin > 1 && (varargin{1}(1)~='-'), load_args = varargin(1); 
  else load_args = {'?'}; 
  end
  
  if any(named('-noise')), load_args = [load_args {'-noise', get_('-n')}]; end    
  if any(named('--l')), load_args = [load_args get_('--l')]; end
  if ~any(named('-no-f')), load_args = [load_args {'-filt'}]; end
  
% if ~exist('data','var'), load_ECAPs('?'); end
% load_ECAPs('-ab','-index',1);
  
  load_ECAPs(load_args{:}); 
  if any(named('-load'))
    assignin('caller','data',data) %#ok<NODEF>
    assignin('caller','list',list)
    assignin('caller','axon_info',axon_info)
    assignin('caller','max_current',max_current)
    return
  end
end

%% Compute example fit (to superimpose on in-vivo data)

% fit_ROI = [0.65 1.2];
% fit_ROI = [3 7]; 
fit_ROI = [4.5 9]; 

if any(named('-roi')), fit_ROI = get_('-roi'); end

roa = data.info.time >= fit_ROI(1) & data.info.time <= fit_ROI(2) ;

for_waves_ = @(f) reshape(arrayfun(@(w) f(w.wave(roa,2)), ...
                                       data.wave(end,:)),[],1);
                                     
vpp = for_waves_( @(y) diff(quantile(y,[0 1])) );      
iua = data.info.stimulus.current; 
% bionics_CL_to_uA =  @(CL) 17.5*100.^(CL/255); % convert, units uA  
bionics_uA_to_CL =  @(uA) 255*log(uA/17.5)/log(100); % convert, units CL
clear for_waves_

CL = round(bionics_uA_to_CL(iua));

if any(named('-ua-r')), iua_range = get_('-ua-r');
  sel = iua >= min(iua_range) & iua <= max(iua_range);
elseif any(named('-CL-r')), CL_range = get_('-cl-r');
  sel = CL >= min(CL_range) & CL <= max(CL_range);
else sel = []; 
end

if ~isempty(sel) % apply CL selection
  
  CL = CL(sel);
  iua = iua(sel);
  vpp = vpp(sel);
  data.wave = data.wave(:,sel);
  for f = fieldnames(data.info)'
    if numel(data.info.(f{1})) ~= numel(sel), continue, end
    data.info.(f{1}) = data.info.(f{1})(sel);
  end
  
end


if ~any(named('-no-sig'))
    
  iua = reshape(iua,[],1);
  vpp = reshape(vpp,[],1);
  
  logit = @(x,p) ((p(2)-p(3))./(1+exp((p(1)-x)./p(4))))+p(3) ; % + p(5)*x;

  u = find(vpp > max(vpp)/2,1); 
  p0 = [ iua(u) max(vpp) min(vpp) iua(u)/4 0];

  LB = [-2*max(iua) 0 0 0 -inf]; 
  UB = @(y) [2*max(iua) 4*max(y) max(y) inf inf];

  fit_opt = optimset('display','off','maxFunEvals',1e5,'maxIter',1e5);

  weight = @(y) abs(y./max(abs(y))).^0.2 ; 
  % weight = @(y) 1; 

  lsq = @(y,p) nanmean((y-logit(iua,p)).^2 .* weight(y) );
  pF = fmincon(@(p)lsq(vpp,p), p0, [],[],[],[],LB,UB(vpp),[],fit_opt);
  pF(5) = 0; 

  clear lsq weight fit_opt LB UB u p0 for_waves_
end


%% waterfall plot for simulated ECAP

if ~any(named('-clf'))
     if isempty(get(gcf,'children'))
          set(gcf,'Position',[100 560 1180 480])
     end
     clf, subplot(1,3,[1 2]), hold on
else cla, hold on
end

view_ROI = [0 12];

if any(named('-view')), view_ROI = get_('-view'); end

roi = data.info.time >= view_ROI(1) & data.info.time <= view_ROI(2) ;
time = data.info.time(roi);

delta_y = 1; 
if any(named('-y')), delta_y = '-y'; end

w_max = max(data.wave(end).wave(roa,2));

[~,~,ci] = ttest(vpp(1:min(3,end)));

threshold = iua(find(vpp > mean([ci(2) max(vpp)]),1));
if threshold == max(iua), threshold= iua(end-1); end
if isempty(threshold), threshold = iua(end-1); end
C = lines(7); 
G = @(v) [v v v]/10; 

color = interp1(log10([min(iua)*0.9; threshold; max(iua)]), ...
                      [(C(6,:)+1)/2; G(1); C(7,:)], log10(iua)); 

rectangle('Position',[fit_ROI(1) -delta_y diff(fit_ROI) ...
                       (1+numel(iua))*delta_y+w_max], ...
                       'FaceColor',[.85 .85 .85],...
                       'EdgeColor','none');

for pp = 1:numel(iua)
  
  plot(time,  data.wave(end,pp).wave(roi,2) + delta_y*(pp-1), 'Color', ...
              color(pp,:), 'LineWidth', 0.8)
            
  if any(named('-ua-l'))
    text(xlim*[-0.02;1.02],delta_y*(pp-1), sprintf('%0.1f µA',iua(pp)),'Color',color(pp,:),'FontSize',6)
  elseif ~any(named('-no-l'))
    text(xlim*[-0.02;1.02],delta_y*(pp-1), sprintf('CL%d',CL(pp)),'Color',color(pp,:),'FontSize',6)
  end
end

ylim([0 max(ylim)]), tools.tidyPlot

[~,s] = fileparts(data.info.path(1:end-1));
title(s)

if any(named('-clf')), return, end

set(gca,'Position',get(gca,'Position') - [0.05 0 0 0])

%% 

% figure(2), clf, hold on

subplot(4,3,[6 9]), cla reset, hold on

if any(named('-ex'))
  plot(iua,vpp,'k-','LineWidth',1.3,'Color',[.8 .8 .8],'clipping','off')
end
plot(iua,vpp,'k-','LineWidth',1.3)

if any(named('-xl')), xlim(get_('-xl'));    
else
  iunit = 10.^floor(log10(iua(end)));
  
  xlim([0.9*min(iua) ceil(iua(end)/iunit)*iunit])
end
set(gca,'XScale','log'), tools.tidyPlot
set(gca,'YScale','log')

xlabel('Stimulus, µA')
ylabel('ECAP V_{PP}, µV')

if ~any(named('-no-sig'))
  x = logspace(log10(iua(1)),log10(iua(end)), 100);
  plot(x, logit(x,pF),'Color',[C(7,:) 0.3],'LineWidth',1.2)
end

%%

if 0  
  %%
  tools.export_H5XP_file(data)
end



function [data, list, axon_info, max_current] = load_ECAPs(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

%% Locate data

if any(named('-dir')), list_path = get_('-dir');
else                   list_path = 'waves~'; 
end

list = dir(tools.file([list_path '\*'])); 
list(cellfun(@(d) all(d=='.'), {list.name})) = []; % remove '..', '.'
list(contains({list.name},'README')) = []; % remove these

b = contains({list.name},'base');
f_ = @(x) [x.folder filesep x.name];

base = list(b);
list = list(~b); 
list(~contains({list.name},'(')) = []; 

%% Select relevent waves and link to baselines 

if any(named('-?')) || any(named('?'))
  
  sel = menu('select ECAP set',{list.name});
  list = list(sel);
  
else  
  sel = ismember({list.name},varargin(cellfun(@ischar,varargin))); 
  if any(sel)
    list = list(sel);
    if ~any(named('-q')), cellfun(@disp, {list.name}), end  
  else 
    % cellfun(@disp, {base.name})
    if ~any(named('-q')), cellfun(@disp, {list.name}), end

    if ~any(named('-all')) && ~any(named('-index'))
      if nargout == 1, data = list; 
      else disp('to load this entire list, load_ECAPs -all')
           disp('to load a selection, load_ECAPs -index [x]')
      end, return
    end
  end
end

clear b

%%  Load necessary data 

axon_info = tools.parse_arguments(varargin,'LOAD','axon'); 

data = struct; % ([]);

if any(named('-index'))
     index = reshape(get_('-index'),1,[]);
else index = 1:numel(list); 
end

for ii = index
  
  if any(named('--d')) % data passed in as an input arg (eg for filter experiment)
    data = get_('--d'); 
    
    for dd = 1:numel(data) % reset wave to w_raw (if exists)
      if ~isfield(data(dd).wave,'w_raw'), continue, end
      for pp = 1:numel(data(dd).wave)                
        data(dd).wave(pp).wave = data(dd).wave(pp).w_raw; 
      end
    end
    
    break, 
  end
  
  %%
  if numel(base) > 1, b = ii; else b = 1; end
  
  mcw_args = {'-chan',[1 2; 3 4],'-quick', '-noise', 0.0275};  
  if any(named('-filter')), mcw_args{end} = 0.275; end
  if any(named('-noise')),  mcw_args{end} = get_('-noise'); end  
  
  if any(named('-type')), mcw_args = [mcw_args {'-type',get_('-type')}];       %#ok<AGROW>
  elseif ~any(named('-class')), mcw_args = [mcw_args {'-type',{[1 2 3 4]}}]; %#ok<AGROW>
  end 
  
  
  if any(named('-args')), more_args = get_('-args'); % generic passthrough
    if iscell(more_args), mcw_args = [mcw_args more_args]; %#ok<AGROW>
    else                  mcw_args = [mcw_args {more_args}]; %#ok<AGROW>
    end
  end

  disp(tools.file('T',f_(list(ii))))
  
  if any(named('-by-f')), rng('default'); end
  
  [data(index == ii).wave, data(index == ii).info] = ...
            models.compose_waves([f_(list(ii)) '\'],mcw_args{:});
          
  if any(named('-by-f'))
    nF = size(axon_info.nerve.coeffs,3); 
    for f_id = 1:nF
      rng('default'); 
     [data(index == ii,f_id+1).wave, data(index == ii,f_id+1).info] = ...
                 models.compose_waves([f_(list(ii)) '\'],mcw_args{:}, ...
                 '-fascicle',f_id);  
    end    
    data(index == ii,:) = data(index == ii,[2:(nF+1) 1]);
  end
          
end

max_current = min(cellfun(@(u) max(u.stimulus.current), {data(:).info})); 

%% Apply filter if requested

if any(named('-filt'))
  data = arrayfun(@(d) apply_filter(d,varargin{:}), data); 
end




%% and return

if any(named('-index')), list = list(index); end

if nargout == 0
  
  assignin('caller','data',data)
  assignin('caller','list',list)
  assignin('caller','axon_info',axon_info)
  assignin('caller','max_current',max_current)
  clear
end

return



function data = apply_filter(data,varargin)
%%
named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

f_opts = [500 3000 4]; % BI default filtering

if any(named('-hz')), f_opts = get_('-hz'); end

if numel(f_opts) == 1, f_opts = [0 f_opts 4];
elseif numel(f_opts) == 2, f_opts(3) = 4;
end

n = abs(f_opts(3)); 
dt = mean(diff(data.info.time)) / 1e3;
  
if f_opts(1) == 0
  fd = fdesign.lowpass('N,F3db',n,f_opts(2),1/dt);
elseif f_opts(2) == 0    
  fd = fdesign.highpass('N,F3db',n,f_opts(1),1/dt);
else     
  fd = fdesign.bandpass('N,F3dB1,F3dB2',n,f_opts(1),f_opts(2),1/dt);    
end
  
fprintf('Filtering wave, %d - %d Hz\n', f_opts(1:2))

Hd = fd.design('butter'); 
set(Hd,'Arithmetic','double');  
[b,a] = sos2tf(Hd.sosMatrix, Hd.ScaleValues);

%%
[data.wave.w_raw] = deal([]); 

for ii = 1:numel(data.wave)
  
  if isempty(data.wave(ii).w_raw)
    data.wave(ii).w_raw = data.wave(ii).wave;
  end
  
  wave = data.wave(ii).w_raw;

  if f_opts(3) > 0, wave = flipud(wave); end  
  
  if any(named('-zerolag'))
    for pp = 1:size(wave,2)  
        wave(:,pp) = cast(filtfilt(b,a,double(wave(:,pp))),'like',wave);  
    end  
  else
    for pp = 1:size(wave,2)  
        wave(:,pp) = cast(filter(b,a,double(wave(:,pp))),'like',wave);  
    end
  end
  if f_opts(3) > 0, wave = flipud(wave); end
  
  if any(~isfinite(wave(:)))    
    [ty,pp] = ind2sub(size(wave),ii);    
    warning('NaN or Inf in wave_stim-%d_axons-%d', pp, ty)
  end
  
  if any(named('-debug-filter'))
    %%
    clf
    plot(data.info.time,data.wave(ii).w_raw(:,end))
    hold on, axis(axis)
    plot(data.info.time,wave(:,end))
  end
  
  data.wave(ii).wave = wave;   
end





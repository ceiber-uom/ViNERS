

function preview_wave(filename,varargin)
% Basic visualisation of waves output by models.nerve_recording and
% models.ECAP_recording
%  
% -f [list] : if wave includes multiple fascicles, filter summation 
%             (Default: sum across all)
% -roi [tMin tMax] : set display ROI, default is ±40 ms.
% -pair [list], -chan [list], -elec [list] : define electrode montages, 
%              default is sequential bipolar pairs.
% -zscale [0.5] : separation between wave baselines, scaled by wave Vpp
% -trend : include linear trend in wave display, otherwise linear trend 
%          is removed using tools.detrend_wave
% 
% V0.1 CDE 30 June 2021

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if nargin == 0, filename = '-prompt'; end
if strcmp(filename,'-prompt') || isnumeric(filename) % output file
  wav_folder = tools.file('get','sub~/waves/','newest');
  
  if isnumeric(filename)
    list = dir([wav_folder '*.mat']);
    filename = [list(filename).folder filesep list(filename).name];
  else
   [filename,wav_folder] = uigetfile('*.mat',[],wav_folder);
    filename = [wav_folder filename];
  end
end

load(filename,'waves','time','options')

time_range = (max(abs(time)) - 25);
bipolar = reshape(1:size(waves,2), 2, [])';
dy_factor = 0.5; 

if any(named('-roi')), time_range = get_('-roi'); end
if any(named('-pair')), bipolar = get_('-pair'); end % vague synonyms
if any(named('-elec')), bipolar = get_('-elec'); end
if any(named('-chan')), bipolar = get_('-chan'); end
if any(named('-z')),    dy_factor = get_('-z'); end
  
if ischar(time_range), time_range = str2double(time_range); end
if numel(time_range) == 1, time_range = [-1 1]*time_range; end

roi = (time >= min(time_range) & time <= max(time_range)); % units of ms

if ~any(named('-tr'))  
  waves = tools.detrend_wave(waves,time,roi);
end

if numel(size(waves)) == 4

  fascicles = (1:size(waves,3));
  if any(named('-f')), fascicles = get_('-f'); end
  waves = squeeze(sum(waves(:,:,fascicles,:),3)); 
end

axon_type = options.class; 
aff = [options.afferent] == 1; 
axon_type(aff)  = cellfun(@(s) sprintf('aff. (%s)',s), axon_type(aff),'unif',0);
axon_type(~aff) = cellfun(@(s) sprintf('eff. (%s)',s), axon_type(~aff),'unif',0);


for type = 1:numel(axon_type)

  
  wb = sum(waves(:,:,type,:),4);
  wb = wb(roi,bipolar(:,1)) - wb(roi,bipolar(:,2));

  dy = 2 * diff(quantile(wb(wb>0), [0.01 0.99]));
  if isnan(dy), dy = 1; end

  subplot(2,2,type), cla, hold on
  C = lines(max(7,size(wb,2)));

  for cc = 1:size(wb,2)
  
    w0 = (cc-1)*dy*dy_factor;
    plot(time(roi), wb(:,cc) + w0,'Color',C(cc,:))
    text(max(time(roi)) + 1, w0, ['E' sprintf('%d',bipolar(cc,:))],...
                                 'Color',C(cc,:),'FontSize',14)  
  end
  
  axis tight, tools.tidyPlot    
  title(axon_type{type})
end

linkaxes(get(gcf,'Children'),'x')
[~,filename] = fileparts(filename);
tools.suptitle(strrep(filename,'_','\_'))
function preview_stimulus(s, varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
% get_ = @(v) varargin{find(named(v))+1};

if nargin < 1, s = ''; end

if ~isstruct(s), 
  if any(s == '~'), s = tools.file(s); end
  if ~exist(s,'file'), s = tools.file('stim~\*.mat','prompt'); end
  s = load(s,'stimulus');    
end

if isfield(s,'stimulus'), s = s.stimulus; end
if isfield(s, 'timing')
  s.t = s.timing;
  s.p = s.pattern;
end

T0 = s.t(1); 

if ~any(named('-delay')), time = s.t - s.t(1);
else time = s.t;
end
time = reshape([time; time],[],1); 

stim = s.p; 
stim = [0*stim(1,:); stim(sort([1:end  1:end-1]),:)];

if any(named('-roi')), roi = get_('-roi'); 
else roi = [0 100] - T0 ; 
end

t_stim = time([1 end]); 

time = [min(roi); time; max(roi)]; 
stim = [0*stim(1,:); stim; 0*stim(1,:)];

if ~any(named('-cla')), cla, end
hold on, plot(time, stim,'LineWidth',1.1)

axis tight, tools.tidyPlot
xlabel('time, ms'), set(gca,'YColor','none')
if any(named('-roi')), xlim(roi)
else xlim(t_stim' * [1.02 -0.02; -0.02 1.02])
end


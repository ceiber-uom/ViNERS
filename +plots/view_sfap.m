function view_sfap(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

% files or gather from input
if nargin > 0, file = varargin{1}; else file = {}; end
if isempty(file) || file(1) == '-',
   file = {'sensitivity*.mat'}; 
   if any(named('-af')), file = file([2 1]); end
  [file,folder] = uigetfile(file,[],tools.file('sfap~\'));
  if all(folder == 0), return, end % cancelled
  disp(file)
  file = strcat(folder,file);
elseif strcmpi(file,'auto')
end
if contains(file,'~'), file = tools.file(file); end

%% List of options parallels plots.view_sensitivity

dat = load(file);
dat.filename = file;

if any(named('-anat')), dat.nerve = get_('-anat'); end
if ~isfield(dat.nerve,'fascicles'), dat.nerve.fascicles = dat.nerve.outline; end
nF = size(dat.nerve.fascicles,3);
nE = size(dat.axons(1).sensitivity(1).profile,2);

electrode_list = 1:nE;    % [1 2 9 10]; 
if any(named('-c')), electrode_list = get_('-c'); end
if any(named('-e')), electrode_list = get_('-e'); end

if any(electrode_list > nE) 
  warning('File x only has Y electrodes')
  electrode_list(electrode_list > nE) = []; 
end

opts.fascicle = 1:nF;
if any(named('-f')), opts.fascicle = get_('-f'); end
opts.electrode = electrode_list; 

opts.do_vdc = any(named('-dc')) || any(named('-vdc'));
opts.do_image = any(named('-im')); 
opts.force_mf = any(named('-mf'));
opts.do_bipolar = any(named('-bi'));

opts.axon_type = 'MRG'; % default axon class
opts.xy_res = 61 * opts.do_image; % native resolution 
opts.pk_quantile = 1;
opts.filter_hz = []; % no filter 

if any(named('-ty')), opts.axon_type = get_('-ty'); end
if any(named('-res')), opts.xy_res = get_('-res'); end
if any(named('-pk-q')), opts.pk_quantile = get_('-pk-q'); end
if any(named('-hz')), opts.filter_hz = get_('-hz'); end

clf, choose_sfap_figure(dat,opts)

%%
return

%% Which figure to display? 
function choose_sfap_figure(dat, opts) 

if nargin < 2, opts = struct; end
if ~isfield(opts,'fascicle'),  opts.fascicle = 1; end
if ~isfield(opts,'electrode'), opts.electrode = 1; end
if ~isfield(opts,'xy'),        opts.xy = []; end

if nargin < 1 && evalin('caller','exist(''dat'',''var'')')
    dat = evalin('caller','dat');    
end

if numel(opts.fascicle) > 1 || ...  % Generate multiple-fascicle figure
   numel(opts.electrode) > 1 || opts.force_mf,
      make_panels_multiFascicle(dat,opts)
else  make_panels_interactive(dat,opts);
end



%%
function make_panels_multiFascicle(dat,opts)

nF = max(opts.fascicle); 
nE = numel(opts.electrode);
nP = size(opts.xy,1);
if nP == 1, nP = opts.xy; end

% Standard color-tools for figures
% C = lines(max(f_id)); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

%%

c_sel = strncmp({dat.axons.class},opts.axon_type,length(opts.axon_type));
if ~any(c_sel)
  c_sel = sprintf('|%s',dat.axons.class);
  error('Invalid axon type. Use -ty ''%s''. Requested: "%s"',c_sel(2:end),opts.axon_type)
elseif sum(c_sel) > 1
  c_sel = find(c_sel,1);
  warning('Multiple options for -ty "%s", using "%s"', opts.axon_type, ...
              dat.axons(c_sel).class)
end

axons = dat.axons(c_sel);


clf
xy = cell(nF,1); 

for ee = 1:nE
      
  sfap_vpp = cell(nF,1);
  sfap_vdc = cell(nF,1);
    
  for f0 = 1:numel(opts.fascicle)
    
    ff = opts.fascicle(f0); 
    fok = axons.axon.fascicle == ff; 
      
    if ~isfield(axons,'axon_SFAP')
      error TODO_compute_axons_SFAP_field_from_component_SFAP
    end
    
    xy_f = dat.nerve.fascicles(:,:,ff); 
    
    if isempty(xy{ff})
      if nP ~= 0 || opts.do_image
        
        xyf = unique(xy_f,'rows');
        xy0 = [median(xyf) min(xyf) max(xyf)];
        xy0(3:4) = max(abs(xy0(1:2)-xy0([3 4; 5 6])),[],1);

        [gx,gy] = meshgrid(xy0(1)-xy0(3)*linspace(-1,1,opts.xy_res), ...
                           xy0(2)-xy0(4)*linspace(-1,1,opts.xy_res));

        if opts.do_image, ok = true(size(gx(:))); 
        else ok = inpolygon(gx(:), gy(:), xy_f(:,1), xy_f(:,2));
        end
                         
        xy{ff} = [gx(ok) gy(ok)]; 
      else
        xy{ff} = axons.axon.xy(fok,:);
      end
    end % isempty(xy{ff})
    
    % Check same units 
    if mean(range(xy_f)) > 100*mean(range(xy{ff}))
        xy{ff} = xy{ff} * 1e3; 
    elseif 100*mean(range(xy_f)) < mean(range(xy{ff}))
        xy{ff} = xy{ff} / 1e3; 
    end
    
    e_pair = opts.electrode(ee);
    if opts.do_bipolar
       e_pair = [e_pair mod(e_pair,size(axons.axon_SFAP,2))+1]; %#ok<AGROW>
         wave = sum(axons.axon_SFAP(:,e_pair,fok) .* [1 -1],2);
    else wave = axons.axon_SFAP(:,e_pair,fok);
    end, wave = permute(wave,[1 3 2]);
    
    wave = apply_filter(wave,axons.time,opts); 
    
    sfap_vpp{ff} = (quantile(wave,opts.pk_quantile,1) - ...
                    quantile(wave,opts.pk_quantile-1,1))';
    sfap_vdc{ff} = (sum(wave,1) / mean(diff(axons.time)))';

    if nP ~= 0 || opts.do_image
      
      field = scatteredInterpolant(axons.axon.xy(fok,:),sfap_vpp{ff},'natural','none');
      sfap_vpp{ff} = field(xy{ff}); 
      if opts.do_vdc, field.Values = sfap_vdc{ff}; 
        sfap_vdc{ff} = field(xy{ff});
      end
    end
    
    %% Plot the SFAP summary data
   
    if opts.do_vdc, subplot(nE,2,2*ee-1),
    else            subplot(ceil(nE/2),2,ee),
    end,            hold on, axis equal
     
    if opts.do_image, sfap_vpp{ff} = reshape(sfap_vpp{ff},size(gx,1),[]);
       imagesc(gx(1,:),gy(:,1),sfap_vpp{ff})
       plot(xy_f(:,1),xy_f(:,2),'-','Color',[0 0 0 0.6])
       colormap(gca,[.9 .9 .9; magma(256)])
    else
       plot(xy_f(:,1),xy_f(:,2),'-','Color',[0 0 0 0.3])
       scatter(xy{ff}(:,1),xy{ff}(:,2),[],sfap_vpp{ff}(:),'.')
       colormap(gca,tools.magma)
    end
     
    if ~opts.do_vdc, continue, end
    subplot(nE,2,2*ee), hold on, axis equal
    
    if opts.do_image, sfap_vdc{ff} = reshape(sfap_vdc{ff},size(gx,1),[]);
      imagesc(gx(1,:),gy(:,1),sfap_vdc{ff})
      plot(xy_f(:,1),xy_f(:,2),'-','Color',[0 0 0 0.6])
      colormap(gca,[.9 .9 .9; magma(256)])
    else
      plot(xy_f(:,1),xy_f(:,2),'-','Color',[0 0 0 0.3])
      scatter(xy{ff}(:,1),xy{ff}(:,2),[],sfap_vdc{ff}(:),'.')
      colormap(gca,'magma')
    end
  end
    
  %% Format axes
  
  if opts.do_vdc, subplot(nE,2,2*ee-1)
  else            subplot(ceil(nE/2),2,ee),
  end
  
  axis tight, ylim([0 max(ylim)]), tools.tidyPlot
  vals = quantile(reshape(cat(1,sfap_vpp{:}),[],1),[.01 .99]);
  if numel(unique(vals)) == 1 || any(isnan(vals))
    vals = [nanmin(cat(1,sfap_vdc{:}))-1e-3 ...
            nanmax(cat(1,sfap_vdc{:}))+1e-3];
  end, caxis(vals)  
  
  ylabel(['\bfE' sprintf('%d',e_pair)]);
  
  if ee == 1, title('V_{PP} (µV)'), end

  if ~opts.do_vdc, pause(0.02), continue, end
  
  ch = colorbar; ch.Position(1) = 0.424;
  
  subplot(nE,2,2*ee), ylim([0 max(ylim)]), tools.tidyPlot
  vals = quantile(cat(1,sfap_vdc{:}(:)),[.01 .99]);
  if numel(unique(vals)) == 1 || any(isnan(vals))
    vals = [nanmin(cat(1,sfap_vdc{:}(:)))-1e-3 ...
            nanmax(cat(1,sfap_vdc{:}(:)))+1e-3];
  end, caxis(vals)
  
  ch = colorbar; ch.Position(1) = 0.904;
  if ee == 1, title('V_{DC} (µV)'), end
  
  pause(0.02);
  
end


tools.suptitle(dat.axons(c_sel).class)
%%
return




%%
function make_panels_interactive(dat,opts)


ff = opts.fascicle(1); 
ee = opts.electrode(1);
if opts.do_bipolar
   ee = [ee mod(ee,size(dat.axons(1).axon_SFAP,2))+1]; 
end

nC = numel(dat.axons); 
sfap.vpp = []; 
sfap.wave = []; 
sfap.xy = []; 

for ty = 1:nC

  subplot(nC,3,3*ty), cla, hold on
  axons = dat.axons(ty); 
  fok = axons.axon.fascicle == ff; 
      
  if ~isfield(axons,'axon_SFAP')
    error TODO_compute_axons_SFAP_field_from_component_SFAP
  end
    
  if opts.do_bipolar
       wave = sum(axons.axon_SFAP(:,ee,fok) .* [1 -1],2);
  else wave = axons.axon_SFAP(:,ee,fok);
  end, wave = permute(wave,[1 3 2]);
  
  wave = apply_filter(wave,axons.time,opts); 
  
  plot(axons.time, median(wave,2),'Color',[.4 .4 .4 .4],'LineWidth',1.2, ...
          'Clipping','on')
        
  [~,sel] = min(sum((wave-median(wave,2)).^2));

  plot(axons.time, wave(:,sel),'Color',dat.pop(ty).color, ...
        'LineWidth',1,'Clipping','on')
    
  d = 20 * axons.axon.diameter(axons.axon.subtype); 
  if numel(d) > 1, d = d(fok); end


  sfap(ty).vpp = (quantile(wave,opts.pk_quantile,1) - ...
                  quantile(wave,opts.pk_quantile-1,1))';
  sfap(ty).wave = wave;
  sfap(ty).xy = axons.axon.xy(fok,:); 
  sfap(ty).diam = d; 

  tools.tidyPlot, title(axons.class)
  
  set(gca,'UserData',ty,'ButtonDownFcn', ...
           @(a,b) update_SFAP_anatomy(a,b,sfap(ty)));
         
end

h = get(gcf,'Children'); 
linkaxes(h,'xy')

subplot(1,3,[1 2]), cla, hold on
xy_f = dat.nerve.fascicles(:,:,ff); 

if mean(range(xy_f)) > 100*mean(range(sfap(ty).xy))
    xy_f = xy_f / 1e3;
elseif 100*mean(range(xy_f)) < mean(range(sfap(ty).xy))
   xy_f = xy_f * 1e3; 
end
    
    
fill(xy_f(:,1),xy_f(:,2),[.8 .8 .8],'EdgeColor',[.4 .4 .4],'LineWidth',1.2)


w = 0.2;
xy = sfap(ty).xy; 
scatter(xy(:,1), xy(:,2), sfap(ty).diam, sfap(ty).vpp, 'o','filled')
colormap(tools.magma), caxis(quantile(sfap(ty).vpp,[0.02 0.98]))
axis image, tools.tidyPlot, ylabel(colorbar,'V_{PP}')
axis([min(xy) max(xy)] * [1+w -w 0 0; 0 0 1+w -w; -w 1+w 0 0; 0 0 -w 1+w])
title(axons.class)

plot(xy(sel,1), xy(sel,2),'s','Color',dat.pop(ty).color, ...
        'MarkerSize',10,'LineWidth',1.3,'UserData',ty+0.1)
set(get(gca,'Children'),'hittest','off')
set(gca,'ButtonDownFcn',@(a,b) update_SFAP_panels(a,b,sfap,dat));


return

%%

% the concept is subplot(1,3,[1 2]) = cross-section of Vpp (or vdc if set)
% subplot(2,3,3) and ,6) are SFAP volts-time traces (update from Left/Right
% click of the main image) 

function update_SFAP_panels(hobj,event,sfap,dat)

% xy = [hobj.Children(2).XData(:) hobj.Children(2).YData(:)];
% [~,sel] = min(sum((event.IntersectionPoint(1:2) - xy).^2,2));
h = flipud(hobj.Parent.Children(2).Children);


for ty = 1:numel(sfap)
  
  ax = findobj(hobj.Parent,'UserData',ty);
  
  xy = sfap(ty).xy;  
  [~,sel] = min(sum((event.IntersectionPoint(1:2) - xy).^2,2));
  
  ax.Children(1).YData = sfap(ty).wave(:,sel);
  
  if ty + 0.1 == h(3).UserData, 
    h(3).XData = xy(sel,1);
    h(3).YData = xy(sel,2);
  end
  
end

return

function update_SFAP_anatomy(hobj,~,this)

h = flipud(hobj.Parent.Children(2).Children);
h(2).CData = this.vpp;
h(2).XData = this.xy(:,1);
h(2).YData = this.xy(:,2);
h(2).SizeData = this.diam; 
h(3).Color = hobj.Children(1).Color;

h(3).UserData = hobj.UserData + 0.1; 

caxis(h(1).Parent, quantile(this.vpp,[0.02 0.98]))

return


%% Signal processing functions

function wave = apply_filter(wave,time,opts)

if ~isfield(opts,'filter_hz'), return, end
persistent b a
if isempty(wave), b = []; a = []; return, end
if isempty(a),
  
  default_order = 2; 
  f_opts = opts.filter_hz;
  if isempty(f_opts), return, end
  if ischar(f_opts), f_opts = [500 3000 default_order]; end % default
  
  if numel(f_opts) == 1, f_opts = [0 f_opts default_order];
  elseif numel(f_opts) == 2, f_opts(3) = default_order;
  end

  n = abs(f_opts(3)); 
  dt = mean(diff(time)) / 1e3;

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
end


for pp = 1:size(wave,2)  
    wave(:,pp) = cast(filtfilt(b,a,double(wave(:,pp))),'like',wave);  
end  

return

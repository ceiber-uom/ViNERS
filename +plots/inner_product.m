
function inner_product(varargin)
% plots.inner_product
% 
% This runs models.nerve_recording for a single spike and produces a plot
% illustrating how the membrane current and sensitivity funciton combine to
% produce a single fibre action potential (SFAP). This is particularly
% helpful for non-straight axon geometries
% 
% See also models.nerve_recording, models.axon_sfap
% 
% Options:
%
% -file 'file.mat'   : specify spatial sensitivity function
% -axons 'axons.mat' : specify axon population for SFAP calculation,
%                      (default: newest axons file in axons~/)
% -currents [folder] : specify root folder for membrane currents.
% -model [model]     : axon model (Gaines, MRG, Sundt, or Tigerholm)
%                       strictly speaking -model is not necessary, but lets 
%                       you extend ViNERS with other axon models. 
% -fascicle [#]      : which fascicle to compute SFAPs and visualise? 
%                       --f [#] is a synonym.
% -elec [list]       : set electrodes (default: sequential bipolar pairs),
%                       list is an [n_electrodes x 1] (monopolar) or 
%                       [n_channels x 2] (bipolar) list of electrodes. 
%                       Synonyms: -pair [], -chan []
% -axon [#]          : which axon model should be used to compute SFAPs? 
%                      (default: use the median axon model).
% 
% More Options: 
% -index [list]   : select indices from the grid to display. By default,
%                    the axon producing the largest SFAP and the axon
%                    closest to the fascicle centre are displayed
% -field {cell array of ScatteredInterpolant} : use a modified set of 
%                                               sensitivity functions
% -nX [n]         : size of the grid (n x n) of axons to simulate
% -delta-xy [x y] : displace axons by specified amount (in µm)
% -distortion [#] : apply a systemic distortion to the recorded velocities,
%                    temporal current width, and/or spatial current width. 
%                   -debug-distortion illustrates the applied distortion.
% -ref [value]    : reference x-axis value for spike-times in raster
% -pad [5]        : change before/after padding for SFAP calculation
% -fs [30 kHz]    : set sampling rate for recorded SFAP (default 30 kHz)
% -time [50 ms]   : set time ROI window for SFAP
% -q, -quiet      : suppress most output to console. 
% -debug-units    : set up the calculation and display axon trajectories,
%                    useful for debugging units issues.
% -smooth-n [5]   : apply N-point smooth to computed SFAP 
%                   (default: 5-point smooth).  
% -smooth-k [kernel] : apply convolutional filter to computed SFAP.
% 
% 
% V0.3 - CDE - updated for modern ViNERS finally

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('-q')), printf = @(varargin) []; else printf = @fprintf; end

%% Revised model set-up

printf('Running plots.%s ... \n', mfilename);

[EM,AX] = tools.parse_arguments(varargin, 'LOAD', ...
                          'eidors','eidors~/sens*.mat', 'axons');

nE = size(EM.Fascicle1.pot,1); 

opts = struct;
opts.axons_folder = AX.folder; 
opts.verbose = ~any(named('-q'));
opts.use_parallel = ~any(named('-no-p'));

nerve = AX.nerve;
pop = AX.pop;

for f = fieldnames(EM.utils)' % Create utility local functions
  EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');
  eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
end

if any(named('-model')),    axon_type = get_('-model');
elseif any(named('MRG')),   axon_type = 'MRG';
elseif any(named('sundt')), axon_type = 'Sundt';
elseif any(named('tiger')), axon_type = 'Tigerholm';
elseif any(named('gai')),   axon_type = 'Gaines';
else axon_type = 'MRG'; % {'Gaines' ; 'Sundt' ; 'MRG' }
end

axon_file_ = @(varargin) [opts.axons_folder axon_type ...
                          filesep sprintf(varargin{:})]; 

f_id = 1; % electrode and fascicle

if any(named('-fa')), f_id = get_('-fa'); % fascicle
elseif any(named('--f')), f_id = get_('--f');
end

%% Get sensitivity functions 

elec_sensor = cell(nE,1); 

x = z_(fac_(f_id)); 
y = y_(fac_(f_id));

if any(named('-field')), 
  elec_sensor = reshape(get_('-field'), [], 1);
  nE = numel(elec_sensor); 
else
  for ee = 1:nE
    elec_sensor{ee} = scatteredInterpolant(x, y, x_(fac_(f_id)), ...
                           EM.(sprintf('Fascicle%d',f_id)).pot(ee,:)', ...
                                             'natural','none');
  end
end

%% Get Distortion if specified

if any(named('-dis'))
  opts.distortion = get_('-dis'); 
  if size(opts.distortion,1) == numel(pop)
    opts.distorion_by_type = opts.distortion; 
  end
  % options.distortion. [ v z t ] is read as a matrix or structure by
  % models.spike_to_wave.
  
  if any(named('-debug-dis')), opts.plot_image = true; end  
end

%% Make simplified axon grid in fascicle

f_outline = nerve.fascicles(:,:,f_id)';

n = 6;
if any(named('-nx')), n = get_('-nx') + 1;
elseif any(named('-ng')), n = ceil(sqrt(get('-ng'))) + 1; 
end

mk_vector_ = @(u) conv(linspace(min(u(:)),max(u(:)),n),[1 1]/2,'valid');

% x = mk_vector_(elec_sensor{1}.Points(:,1));
% y = mk_vector_(elec_sensor{1}.Points(:,2));

x = mk_vector_(f_outline(1,:));
y = mk_vector_(f_outline(2,:));

[y,x] = meshgrid(y,x);

clear n

% ok = ~isnan(elec_sensor{1}(x(:),y(:),0*y(:)));
if any(named('-debug-units'))
  %%
  
  axons.xyz = permute(tools.from_trajectory(EM,nerve,[x(:) y(:)]),[1 3 2]);

  % axons.xyz = permute(tools.from_trajectory(EM,nerve,nerve.fascicles),[1 3 2]);
  
  clf
  plot3(elec_sensor{1}.Points(:,3),elec_sensor{1}.Points(:,2), ... 
        elec_sensor{1}.Points(:,1),'.')
  hold on

  plot3(axons.xyz(:,1,1),axons.xyz(:,1,2),axons.xyz(:,1,3),'o','LineWidth',1.2)
  view(90,0), C = lines(7);
  
  for ii = 1:size(axons.xyz,1)
      plot3(axons.xyz(ii,:,1),axons.xyz(ii,:,2),axons.xyz(ii,:,3),'-','Color',[C(2,:) 0.7])
  end
  
  
  if any(named('-and-plot')), figure;
  else return
  end  
end

axons.xyz = tools.from_trajectory(EM,nerve,[x(:) y(:)],'-f',f_id);

ok = elec_sensor{1}(axons.xyz(:,3,:),axons.xyz(:,2,:),axons.xyz(:,1,:));
ok = mean(isnan(ok),3) < 0.2;

axons.xy_grid = [x(ok) y(ok)];
axons.xyz = axons.xyz(ok,:,:);

clear rotate c_xy f ee fn fp d gg x y ok

%% Compute intracellular-to-extracellular relationship

% propegation is in +Z direction 

fs = 30; % 24.414; % kHz 
if any(named('-sample-rate')), fs = get_('-sample-rate'); 
elseif any(named('-fs')),      fs = get_('-fs'); 
end

time_span = 50;
if any(named('-time')),time_span = get_('-time'); end

time = 0:(1/fs):(time_span);
time = [-fliplr(time) time(2:end)]; % ms

% opts.jitter = 0;
% opts.fraction = 1; 
opts.reference_z = 0; % median(EM.model.nodes(EM.model.electrode(3).nodes,1)); 
opts.class = axon_type;
opts.padding = 5;

if any(named('-pad')), 
  opts.padding = get_('-pad');
end

if any(named('-ref')) z_ref = get_('-ref'); end
if ischar(z_ref), z_ref = str2double(regexp(z_ref,'\d+','match','once'));
    z_ref = median(EM.model.nodes(EM.model.electrode(z_ref).nodes,1)); 
end


opts.get_all = true;

sel = [pop.myelinated];
if ismember(axon_type,{'Sundt','Tigerholm'}), sel = ~sel; end
axons.size_sample = cat(1,pop(sel).size_sample);
axons.size_sample = cat(1,pop(sel).size_sample);
axons.fibre_diam  = cat(1,pop(sel).fibre_diam);
axons.g_ratio     = cat(1,pop(sel).g_ratio); 

mafd = median(axons.fibre_diam); 

if isempty(axons.g_ratio)
     [~,sel] = min(abs(axons.fibre_diam - mafd)); 
else [~,sel] = min((axons.fibre_diam - mafd).^2 + mafd * ...
                   (axons.g_ratio - median(axons.g_ratio)).^2 ); 
end
gg = axons.size_sample(sel); 
if any(named('-ax')), gg = get_('-ax'); end

printf('Running models.spike_to_wave ... \n');
[~,V] = models.spike_to_wave(gg,time,elec_sensor,axons.xyz,opts);
% beep % This can take a minute, beep to get my focus back

%% Figure illustrating inner product

% Standard color-tools for figures
C = lines(max(7,nE)); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

% nXY = size(V,3); 
% example_axon = [2 nXY-1]; % [2 size(axonIndex.xy_grid,1)];

if any(named('-smooth-k')), 
     smoothing_kernel = reshape(get_('-smooth-k'),[],1);
elseif any(named('-smooth-n')), 
     smoothing_kernel = ones(get_('-smooth-n'),1);
else smoothing_kernel = ones(5,1); 
end

e_id = 1; % only used for selection of example[1]

% this code will however also work for bipolar e_id
sfap_vpp = convn(permute(V(:,e_id,:),[1 3 2]),smoothing_kernel,'same');
if size(sfap_vpp, 3) > 1
  sfap_vpp = sfap_vpp(:,:,1) - mean(sfap_vpp(:,:,2:end)); 
end
sfap_vpp = range(sfap_vpp);

[~,example_axon(1)] = max(sfap_vpp);
[~,example_axon(2)] = min(sum((axons.xy_grid-mean(axons.xy_grid)).^2,2));

example_chan = reshape(1:nE,[],2)';

if any(named('-id')),      example_axon = get_('-id');
elseif any(named('-ind')), example_axon = get_('-ind');
end

if any(named('-chan')),     example_chan = get_('-chan');
elseif any(named('-elec')), example_chan = get_('-elec');
elseif any(named('-pair')), example_chan = get_('-pair');
end


%% Render figure

clf
subplot(1,5,1), hold on % upper left panel (anatomy graphic)
plot(f_outline(1,:),f_outline(2,:),'-','Color',G(5),'LineWidth',1.4)
rectangle('Position',[min(xlim) -0.01 range(xlim) 0.01 ],'FaceColor',G(3),'EdgeColor','none')

axis equal tight, tools.tidyPlot, set(gca,'XColor','none')

h = gobjects(0); 

for row_id = 1:numel(example_axon), id = example_axon(row_id);     
    h(row_id) = text(axons.xy_grid(id,1)+0.006, ...
                  axons.xy_grid(id,2),char('A'+row_id-1), ...
                  'FontSize',14, 'Color',G(3));
end

plot(axons.xy_grid(:,1),axons.xy_grid(:,2),'o','Color',G(4), ...
            'MarkerFaceColor',G(8),'MarkerSize',6)
set(gca,'Position',get(gca,'Position') - [0.08 0 0 0])




% selected example - from V computation 
d = load(axon_file_('n%03d_out.mat',gg));
[~,t0] = max(d.V_example(:,2));
dz = mean(diff(d.length)); 
mem_voltage = @(z) interp1(d.time-d.time(t0),  ...
                           d.V_example(:,2), time-z*d.dt_dx/dz, ...
                            'Linear',d.V_example(1));
mem_current = @(z) interp1(d.time-d.time(t0),  ...
                           d.I_membrane(:,2), z*d.dt_dx/dz-time, ...
                            'Linear',d.I_membrane(1));

nEX = numel(example_axon);      

for row_id = 1:nEX, id = example_axon(row_id);     
    %% Plot example Sensitivity cirve and I_membrane
    h(end+1) = subplot(nEX+1,5,5*row_id-[3 2]);
    cla, hold on
    
    % z = linspace(-5,5,201);
    % xy = axons.xy_grid(id,:);     
    % xyz = [xy 0*xy(:,2)] + z' * [0 0 1];    
    
    xyz = permute(axons.xyz(id,[3 2 1],:),[3 2 1]);
    z = xyz(:,3)';
    if numel(z) < 30,       
      xyz = interp1(z',xyz,linspace(z(1),z(end),201)'); 
      z = xyz(:,3)';
    end
    
    yl = []; 
    
    mc_img = cell(size(elec_sensor));
    vm_img = cell(size(elec_sensor));

    for ee = 1:size(example_chan,1)
      
        pair = example_chan(ee,:);
        
        if numel(pair) > 1, pot = elec_sensor{pair(1)}(xyz) - ...
                                  elec_sensor{pair(2)}(xyz); 
        else                pot = elec_sensor{pair(1)}(xyz);
        end
        ok = ~isnan(pot); % Patch over discontinuities
        pot = interp1(z(ok),pot(ok),z,'linear','extrap');
        
        % dz = mean(diff(d.length))
        
        fill(z([1 1:end end]),[0 pot 0],C(ee,:), ... 
                    'EdgeColor','none','FaceAlpha',0.1,'UserData',ee)
        yl = [min([pot yl]) max([pot yl])];
        [~,pmax(ee)] = max(abs(pot));         %#ok<AGROW>
        
    end
    ylim(yl*[1.04 0; -0.04 1])
    
    for ee = 1:size(example_chan,1)
      
        pair = example_chan(ee,:);
        eidx = cellfun(@(n) reshape(n,1,[]), ...
                        {EM.model.electrode(pair).nodes}, 'unif',0);
        eidx = [eidx{:}];
        z_ch = nanmedian(EM.model.nodes(eidx,1));
        z_eq = time/d.dt_dx * mean(diff(d.length));
        vm = mem_current(z_ch);
        
        ok = (z_eq > -5 & z_eq < 5);
        vm = (vm-min(vm))/max(range(vm),1e-3) * range(yl) + yl(1);
        plot(z_eq(ok),vm(ok),'-','Color',C(ee,:),'UserData',[ee id yl])
    end
    
    if row_id == 2, xlabel('Distance in space, mm')
    else set(gca,'XTickLabel','')
    end
    ylabel('Sensitivity, µV/µA')
    tools.tidyPlot
        
    set(gca,'Position',get(gca,'Position') - [0.04 0 -0.04 0])    
    xl = xlim;
    
    %% Plot SFAP 
        
    h(end+1) = subplot(nEX+1,5,5*row_id-[1 0]); %#ok<AGROW>
    set(gca,'Position',get(gca,'Position') + [0 0 0.04 0])
    cla, hold on
    
    vm = mem_voltage(0) / 1e3; 
    % vm = (vm-min(vm))/range(vm)*range(yl)/2 + yl*[-0.1;1.1];
    plot(time,vm,'Color',G(6.5),'LineWidth',1.2)
    
    for ee = 1:size(example_chan,1)      
        pair = example_chan(ee,:);
        
        if numel(pair) > 1, meas = V(:,pair(1),id)-V(:,pair(2),id);
        else                meas = V(:,pair(1),id);
        end
        
        plot(time,meas,'Color',C(ee,:),'LineWidth',1.2)
    end
    tools.tidyPlot, axis tight,
    if ~strcmp(axon_type,'Sundt'),xlim([-1.2 1.2]), end
    axis(axis*[1 -0.02 0 0; 0 1.02 0 0; 0 0 1.04 0; 0 0 -0.04 1]);
    
    if row_id == 2, xlabel('time, ms')
    else set(gca,'XTickLabel','')
    end
    set(gca,'YAxisLocation','right')    
    ylabel('signal, µV')

end

% 
h(end+1) = subplot(nEX+1,5,(nEX+1)*5-[3 2]);
set(gca,'Position',get(gca,'Position') - [0.04 0 -0.04 0])

cla, hold on

sel = EM.model.object_id{strcmp(EM.model.object_name,'PDMS')}; 
idx = EM.model.elems(sel,:);
idx = idx(sum(reshape(y_(idx),[],4) >= -10*eps,2) > 2,:);

patch('Faces',idx, 'Vertices',EM.model.nodes(:,[1 3]), ...
      'EdgeColor',G(7),'FaceColor','w','EdgeAlpha',0.2)

a_outline = EM.model.nodes(unique(idx),[1 3]);
idx = convhull(a_outline); 
plot(a_outline(idx,1),a_outline(idx,2),'-','Color', G(5),'LineWidth',1.2)

axis equal tight, xlim(xl), axis off

for ee = 1:nE

    e_n = EM.model.electrode(ee).nodes;
    e_xy = [min(x_(e_n)) min(z_(e_n)) range(x_(e_n)) range(z_(e_n))];    
    rectangle('Position',e_xy,'FaceColor',W(ee,2),'EdgeColor',W(ee,0.6),'LineWidth',1.2)
end

for row_id = 1:2, id = example_axon(row_id);     
    
    % z = linspace(-5,5,201);
    % xy = axons.xy_grid(id,:);     
    % xyz = [xy 0*xy(:,2)] + z' * [0 0 1];
    
    xyz = permute(axons.xyz(id,[3 2 1],:), [3 2 1]);
    
    tex = {char('A'+row_id-1), 'FontSize',14, 'Color',G(3),'Vert','top'};
    if row_id == 2, tex{end} = 'bottom'; end
    
    plot(xyz(:,3),xyz(:,1),'-','Color',G(3),'LineWidth',1.1)
    text(xyz(end,3),xyz(end,1),tex{:});
end

clear xyz xy z sel idx pot ee e_n e_xy ans vm eidx

suptitle(strrep(EM.filename,'_','\_'));



%% Another snippet of plot code ? 


% in < for ee = example_chan > 
% 
%         ok = abs(time) < mean(abs(time));         
%         [~,tmax(ee)] = max(abs(V(:,ee,id)) .* ok');  %#ok<SAGROW>
%         
%         mc_img{ee} = []; 
%         vm_img{ee} = []; 
%         for pp = 1:numel(z)-1
%             i2v_seg = pot( (pp-1)*nZ+(1:nZ) );            
%             c = current(pp-p0, -14 );
%             mc_img{ee} = [mc_img{ee} c(1,:)];
%             vm_img{ee} = [vm_img{ee} c * i2v_seg];
%         end


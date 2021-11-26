
function list = compare_sensitivity(varargin)
% function [files] = compare_sensitivity([files], ...)
% 
% Compare sensitivity functions in fascicle centre across different meshes.
% If only one file is specified, the comparison is instead across fascicle
%   centres within a mesh. 
% 
% Options
% -fascicle [id] : use specified fascicle, default F#1
% -AF            : compute activating function and show that (for stimulus)
% -bipolar       : compute bipolar sensitivity montage
% -chan [elec]   : show listed electrodes only (Default: all electrodes)
% -elec [elec]   : synonym for -chan
% -no-mesh       : do not show meshes on left hand side of figure (faster)
% -multiplot     : also plot each individual file
%                  (ala plots.view_sensitivity)
% 
% V0.2 CDE 

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

% files or gather from input
if nargin > 0, list = varargin{1}; else list = {}; end
if isempty(list)
  [list,folder] = uigetfile('sensitivity*.mat',[], ...
                              tools.file('sub~\eidors\'), ... 
                             'Multiselect','on');
  if all(folder == 0), return, end % cancelled
  list = strcat(folder,list);
  
end
if isstruct(list) && isfield(list,'name') && isfield(list,'folder'), ...
    p_ = @(x,varargin) [x.folder filesep x.name varargin{:}];
    list = arrayfun(p_,list,'unif',0);
end
if ~iscell(list), list = {list}; end
list(contains(list,'~')) = cellfun(@tools.file, list(contains(list,'~')), 'Unif',0);

tools.setupEIDORS;

nMaps = numel(list);

%% Assuming all files have same electrodes ... 

C = 1-summer(max(7,nMaps)); 
if any(named('-colormap')), C = get_('-colormap'); end

figure(1), clf, lbl = cell(size(list));

fac_id = 1;
if any(named('-f')), fac_id = get_('-f'); end
do_fascicle_comparison = (numel(list) == 1); 
  
if do_fascicle_comparison
  
  fprintf('Loading %s\n',tools.file('T',list{1}))
  EM = load(list{1});
  nMaps = sum(strncmpi(EM.model.object_name,'Fascicle',6)); 
  C = lines(min(7,nMaps));
end

for ff = 1:nMaps
  
  if ff == 1 || ~do_fascicle_comparison

    if ~do_fascicle_comparison
      fprintf('Loading %s\n',strrep(list{ff},tools.file,'~'))  
      EM = load(list{ff});
    end
    
    lbl{ff} = sprintf('(%dx%dx%d mm)',round(range(EM.model(1).nodes(:,[3 1 2]))));
  
    FLAG_stimulus = isfield(EM,'v_extracellular');
    if FLAG_stimulus
      % if FLAG_verbose, 
        if any(named('-AF'))
             disp('Formatting V_e into activating function in fascicles (-AF)')
        else disp('Formatting V_e into fascicles (-STIM)')
        end
      % end
      EM = convert_stim2fascicles(EM);

    elseif any(named('-bi'))
      % if FLAG_verbose, 
        disp('Formatting response into bipolar measurements')
      % end
      EM = convert_mono2bipolar(EM); 
    end

    do_multiMesh = numel(EM.model) > 1; 

    mk_shortcuts(EM,do_multiMesh)

    FascicleX = sprintf('Fascicle%d',fac_id);
    
  else
    
    FascicleX = sprintf('Fascicle%d',ff);
    fac_id = ff; 
  end
  
  nE = size(EM.(FascicleX).pot,1); 
  
  electrode_list = 1:nE; % [1 2 7 8]; 

  if any(named('-chan')), electrode_list = get_('-chan'); end
  if any(named('-elec')), electrode_list = get_('-elec'); end
  
  if any(electrode_list > nE) 
    warning('ViNERS:nElectrodes','%s only has %d electrodes',list{ff},nE)
    electrode_list(electrode_list > nE) = []; 
  end
  
  elec_sensor = cell(nE,1); 

  for ee = 1:nE
    if exist('electrode_list','var') && ~ismember(ee,electrode_list)
      continue
    end
    elec_sensor{ee} = scatteredInterpolant(z_((fac_id)), y_((fac_id)), x_((fac_id)), ...
                                           EM.(FascicleX).pot(ee,:)','natural','none'); 
    
  end

  % Standard color-tools for figures
  W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10;
  
  if isfield(EM.info,'FascicleTrajectory')
    
      alist = dir(tools.file('sub~\axons\*.mat'));    
      load(tools.INPUT_file(alist,list{ff}),'nerve');
      
      xy = mean(nerve.fascicles); 
      xyz = tools.from_trajectory(EM,nerve,xy);
      xyz = permute(xyz(:,[3 2 1],:),[3 2 1]);
      z = cumsum(sqrt(sum(diff(xyz).^2,2))); % length of axon
      z(end+1) = [2 -1]*z([end end-1]); % extend by 1 (because diff)
      z = z - z(find(xyz(:,3) == min(xyz(:,3)),1));
  else   
      xy = [mean(z_((fac_id)))   mean(y_((fac_id)))  0];  
      z = linspace(-8,8,401);
      xyz = xy + z' * [0 0 1];
  end

  for ee = 1:nE

      if exist('electrode_list','var') 
        if ~ismember(ee,electrode_list), continue
        else subplot(numel(electrode_list),1,find(electrode_list == ee))
        end
      else subplot(nE,1,ee)
      end
    
      pot = elec_sensor{ee}(xyz);
      
      % pot = pot ./ max(pot(:));
      
      plot(z,pot,'Linewidth',1.2,'Color',C(ff,:), ... 
                 'userData',sprintf('Elec%d',ee)), hold on
  end
end

h = get(gcf,'Children');

for ii = 1:length(h)
    
    x0 = min(xlim(h(ii))); tools.tidyPlot(h(ii))
    if ii > 1, set(h(ii),'XTickLabel',''), end
    set(h(ii),'TickLength',[1 1]/100)
    text(x0,max(ylim),h(ii).Children(1).UserData, ... 
          'Color',get(gca,'XColor'),'FontSize',14,'PArent',h(ii))
end

set(h,'YLim',[min([h.YLim]) max([h.YLim])])

% lh = legend(lbl);
% legend(cellfun(@(f)regexp(f,'\([^\)]+\)','match','once'),files,'unif',0))
arrayfun(@(n) grid(n,'on'),h)

%% Add mesh on left side
if ~any(named('-nomesh')) && ~any(named('-no-m')) 
  %%
  if isempty(get(gcf,'Children'))  
       subplot_ = @(f) subplot(ceil(numel(list)/2),2,f);    
       % subplot_ = @(f) subplot(numel(files),1,f);    
  else subplot_ = @(f) subplot(numel(list),3,3*f-2);  
    for ee = 1:length(h)
      h(ee).Position([1 3]) = [0.35 0.6]; 
    end
  end

  for ff = 1:numel(list)
    EM = load(list{ff});  
    subplot_(ff)
    
    if numel(EM.model) == 1, show_fem(EM.model)
    else    show_fem(EM.model(fac_id))
    end
    u = get(gca,'Children');
    u(end).EdgeAlpha = 0.1;
    set(gca,'Position',get(gca,'Position') - [5 -1 0 5]/100)  
    title(regexp(list{ff},'\([^\)]+\)','match','once'),'Color',C(ff,:))  
  end

  if numel(get(gcf,'Children')) == nMaps
    h = get(gcf,'children');
    lim = [-1 1]*max(abs([h.XLim])); 
    set(h,'XLim',lim,'YLim',lim,'ZLim',lim)
  end

  % Clean up ticks and add some formatting 

  h = get(gcf,'Children');
  h(arrayfun(@(s) strcmp(s.Type,'legend'), h)) = [];
  h(cellfun(@(s) isempty(s.String),{h.Title})) = []; 
  
  set(h,'XTick',-30:30,'YTick',-30:30,'ZTick',-30:30)
  set(h,'XTickLabel','','YTickLabel','','ZTickLabel','')

  for ii = 1:numel(h)
    h(ii).Position = h(ii).Position+[0 -0.02 0 0.04];
  end
end

%%
if ~any(named('-multip'))
  if nargout == 0, clear, end
  return
end
%% Plot each file (duplicate of fig.1,  plots.view_sensitivity)

figure(2), clf

for ff = 1:numel(list)

  EM = load(list{ff});

  for f = fieldnames(EM.utils)' % Create utility local functions
    EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
    eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
  end
  
  do_multiMesh = numel(EM.model) > 1; 
  
  if do_multiMesh, 
    x_ = @(i) EM.model(fac_id).nodes(i,1);
    y_ = @(i) EM.model(fac_id).nodes(i,2);
    z_ = @(i) EM.model(fac_id).nodes(i,3);
  end
  
  nE = sum(~cellfun(@isempty,{EM.model(1).electrode.name})); 
  elec_sensor = cell(nE,1); 

  
  
  for ee = 1:nE
      elec_sensor{ee} = scatteredInterpolant(z_((fac_id)), ...
                           y_((fac_id)), x_((fac_id)), ...
                               EM.(FascicleX).pot(ee,:)','natural','none'); 
  end

  % Standard color-tools for figures
  C = lines(max(7,nE)); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10;
  
  xy = [mean(z_((fac_id)))   mean(y_((fac_id)))  0];  

  if do_multiMesh, m = fac_id; else m = 1; end
  
  subplot(numel(list),1,ff), hold on

  sel = strcmp(EM.model(m).object_name,'PDMS'); 
  sel = EM.model(m).object_id{sel}; 
  idx = EM.model(m).elems(sel,:);
  idx = idx(sum(reshape(y_(idx),[],4) >= -10*eps,2) > 2,:);
  
  y_gscale = [1 .25 .35];
  if min(EM.(FascicleX).pot(:)) < 0.5  
    y_gscale(3) = max(EM.(FascicleX).pot(:)) + 0.025 - min(z_(idx))*y_gscale(2);
  end

 
  patch('Faces',idx, 'Vertices',EM.model(m).nodes(:,[1 3]) .* y_gscale(1:2) + [0 y_gscale(3)], ...
        'EdgeColor',[.7 .7 .7], 'FaceColor','w','EdgeAlpha',0.2)

  outline = EM.model(m).nodes(unique(idx),[1 3]);
  idx = convhull(outline); 

  plot(outline(idx,1),outline(idx,2)*y_gscale(2) + y_gscale(3), '-', ...
      'Color', [.5 .5 .5],'LineWidth',1.2)


  if isfield(EM.info,'FascicleTrajectory')
    
      alist = dir(tools.file('sub~\axons\*.mat'));    
      load(tools.INPUT_file(alist,list{ff}),'nerve');
      
      xy = mean(nerve.fascicles); 
      xyz = tools.from_trajectory(EM,nerve,xy);
      xyz = permute(xyz(:,[3 2 1],:),[3 2 1]);
      z = cumsum(sqrt(sum(diff(xyz).^2,2))); % length of axon
      z(end+1) = [2 -1]*z([end end-1]); % extend by 1 (because diff)
      z = z - z(find(xyz(:,3) == min(xyz(:,3)),1));
  else   
      xy = [mean(z_((fac_id)))   mean(y_((fac_id)))  0];  
      z = linspace(-8,8,401);
      xyz = xy + z' * [0 0 1];
  end    
    
  for ee = 1:nE

      pot = elec_sensor{ee}(xyz);

      plot(z,pot,'Linewidth',1.2,'Color',C(ee,:)), hold on

      e_n = EM.model(m).electrode(ee).nodes;

      e_xy = [min(x_(e_n)) min(z_(e_n))*y_gscale(2) + y_gscale(3) ... 
              range(x_(e_n)) range(z_(e_n))*y_gscale(2)]; 

      rectangle('Position',e_xy,'FaceColor',W(ee,1), ... 
                'EdgeColor',C(ee,:),'LineWidth',1.2)
    %           ,'Curvature',[1 1]
  end
  


  ok = ~isnan(pot); 
  plot(z(ok),xyz(ok,1)*y_gscale(2) + y_gscale(3),'--','Color',G(3),'LineWidth',1.2)
  set(gca,'DataAspectRatio',y_gscale./min(y_gscale))

  axis tight, tools.tidyPlot
  title(regexp(list{ff},'\([^\)]+\)','match'))
  pause(0.05)
end

h = flipud(get(gcf,'Children'));
set(h,'XLim',[-1.02 1]*max(abs([h.XLim])))
arrayfun(@(n) grid(n,'on'),h)

if nargout == 0, clear, end



%% Utility Functions
function EM = convert_stim2fascicles(EM)
  
  FLAG_multimesh = (numel(EM.model) > 1); 
  named = evalin('caller','named');
  

  if FLAG_multimesh, EM.model = [EM.model{:}];
       nF = numel(EM.model); 
  else nF = sum(contains(EM.model.object_name,'Fasc')) - ...
            sum(contains(EM.model.object_name,'P_Fasc'));
  end
  
  fasc_ = @(n) sprintf('Fascicle%d',n);
  
  % EM.v_extracellular is a list, (nodes x electrodes) of the sequential
  % bipolar stimulus-induced extracellular potential, for ALL nodes 
  % (not just fascicle nodes). 

  % a stimulus file is broken up by fascicle, which each fascicle
  % containing a NODE index and a list of per-node values 
  
  % EM.model.object_id is a list of ELEMENT ids which make up each object
  % in the EIDORS model. Gather each NODE in that set of elements, then use
  % that to emulate the EM.FascicleN.pot/idx 
    
  for ff = 1:nF
    
    if FLAG_multimesh, m = ff; else m = 1; end
    
    sel = strcmpi(EM.model(m).object_name,fasc_(ff));
    if sum(sel) == 0 && FLAG_multimesh, 
      sel = strcmpi(EM.model(m).object_name,fasc_(1));      
    end
    if sum(sel) ~= 1, error('Object %s not found in EM.model.object_name',fasc_(ff)); end
    idx = unique(EM.model(m).elems(EM.model(m).object_id{sel},:));
    
    EM.(fasc_(ff)).idx = idx'; 
    if FLAG_multimesh
         EM.(fasc_(ff)).pot = EM.v_extracellular{ff}(idx,:)';
    else EM.(fasc_(ff)).pot = EM.v_extracellular(idx,:)';
    end
  end
  
  if any(named('-AF'))
    mk_shortcuts(EM,FLAG_multimesh)
    for ff = 1:nF % ACTIVATING FUNCTION (if requested) 
      dx = 1e-3; 
      Ve = scatteredInterpolant(z_(ff), y_(ff), x_(ff), ...
                                EM.(fasc_(ff)).pot(1,:)','natural','none'); 
      for ee = 1:nE

        Ve.Values = EM.(fasc_(ff)).pot(ee,:)';
        xyz = [z_(ff) y_(ff) x_(ff)-dx];  vx0 = Ve(xyz); % have to 2-step
        xyz = [z_(ff) y_(ff) x_(ff)+dx];  vx1 = Ve(xyz);

        EM.(fasc_(ff)).pot(ee,:) = (vx0+vx1-2*Ve.Values)/(dx^2);
      end
    end
  end
  
return
function EM = convert_mono2bipolar(EM)

  FLAG_multimesh = (numel(EM.model) > 1);  

  if FLAG_multimesh, EM.model = [EM.model{:}];
       nF = numel(EM.model); 
  else nF = sum(contains(EM.model.object_name,'Fasc')) - ...
            sum(contains(EM.model.object_name,'P_Fasc'));
  end
  
  nE = size(EM.Fascicle1.pot,1);
  fasc_ = @(n) sprintf('Fascicle%d',n);
  ret_ = @(e) mod(e,nE)+1; % sequentual bipolar pattern

  for ff = 1:nF, mono = EM.(fasc_(ff)).pot;
    for ee = 1:nE
      EM.(fasc_(ff)).pot(ee,:) = mono(ee,:) - mono(ret_(ee),:); 
    end
  end
  for ee = 1:nE % set "stim" structure 
    
    this = struct;
    this.stim_pattern = zeros(nE,1);
    this.stim_pattern(ee) = 1;
    this.stim_pattern(ret_(ee)) = -1;   
    if isempty(EM.model.stimulation)
         EM.model.stimulation = this; 
    else EM.model.stimulation(ee) = this;
    end
    
  end
  
  EM.v_extracellular = false; 
  
return

function mk_shortcuts(EM,is_multimesh)

% if evalin('caller','exist(''x_'',''var'')'), return, end
if nargin == 1, is_multimesh = (numel(EM.model) > 1); end

if is_multimesh
  
  mdl_ = @(n) EM.model(n);
  fac_ =  @(n) EM.(sprintf('Fascicle%d',n)).idx;
  x_ = @(n) EM.model(n).nodes(fac_(n),1); 
  y_ = @(n) EM.model(n).nodes(fac_(n),2); 
  z_ = @(n) EM.model(n).nodes(fac_(n),3);
else  
  mdl_ = @(n) EM.model(1);
  fac_ =  @(n) EM.(sprintf('Fascicle%d',n)).idx; 
  x_ = @(n) EM.model.nodes(fac_(n),1); 
  y_ = @(n) EM.model.nodes(fac_(n),2); 
  z_ = @(n) EM.model.nodes(fac_(n),3); 
    
  % for f = fieldnames(EM.utils)' % Create utility local functions
  %   EM.utils.(f{1}) = strrep(EM.utils.(f{1}),'out','EM');     
  %   eval(sprintf('%s = %s;',f{1},EM.utils.(f{1})));
  % end
end

nE = size(EM.Fascicle1.pot,1);
nF = sum(contains(fieldnames(EM),'Fascicle'));

% Standard color-tools for figures
C = lines(max(7,nE)); W = @(i,v) (C(i,:)+v)/(1+v); G = @(v) [v v v]/10; 

for k = {'mdl_','fac_','x_','y_','z_','nE','nF','C','W','G'}
  assignin('caller',k{1},eval(k{1}));
end


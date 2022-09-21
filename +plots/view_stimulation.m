
function view_stimulation(varargin)
% view_stimulation(filename, ...)
% make 2d cross-sections of stimulation or recording sensitivity profiles
% (a different viewer than plots.view_sensitivity, focus on section planes)
% this also shows the mesh at the cut plane (at the transverse plane) for 
%  additional context
% 
% plots.view_stimulation('-newest', ...) plots most recent file
% plots.view_stimulation('-pdf', ...) makes a PDF from all the files in
%                                         this subject (or -pdf-list, ...)
% 
% Optional arguments: -z [1.875 mm] - location of transverse section plane
%  [default values]   -z auto uses the maximum extracellular potential z
%                     -y [0 mm] - location of longitudinal section plane
%                                 (measured from fascicle centers)
%                     -roi [0.2 mm] - size of transverse view
%                     -elec [1] - stimulation or recording electrode pair
% 
% CDE 24-Nov-2020 v0.3

warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};

if any(named('-pdf'))  
  make_this_PDF(varargin{:}) % Run this as a script
  return
end

if nargin > 0, filename = varargin{1}; 
  if any(named('-ne'))
     filename = tools.file('get','sub~/eidors/stimulus*.mat','newest');
  end  
else filename = tools.file('sub~/eidors/');
  if length(dir([filename 'stim*.mat'])) == 1
    filename = dir([filename 'stim*.mat']); 
    filename = [filename.folder filesep filename.name];
  end
end

if ~isfile(filename) % or exist(filename,'file') == 2  
  
  fp = fileparts(filename); 
  if isempty(fp), fp = tools.file('eidors~/'); end
  [fn,fp] = uigetfile('stimulus*.mat',[],fp);
  if all(fp == 0), clear, return, end
  filename = [fp fn];
  clear fn fp
end

elec_id = 1; % [1 2 7 8]; 
if any(named('-chan')), elec_id = get_('-chan'); end
if any(named('-elec')), elec_id = get_('-elec'); end


fprintf('Viewing %s\nElec %d\n',strrep(filename,tools.file,'~'), elec_id)
load(filename) %#ok<LOAD>


if ~exist('model','var') && exist('m','var'), model = m; clear m, end
multi_geometry_model = iscell(model);
if multi_geometry_model, 
  warning('fascicle_potential_section:multiGeometry', ...
          'Multiple model geometries not fully supported')    
  model = model{1};
end

nE = sum(strncmpi({model.electrode.name},'Elec',4)); 

%% Set up utilities and fascicle list 

if isempty(get(0,'Children')), figure(1), clf
  set(gcf,'Position',[50 200 1380 420])
end

cla reset

C = lines(7); G = @(v) [v v v]/10; 

obj_tet_ = @(n) cat(1,model.object_id{strncmpi(model.object_name,n,length(n))});
obj_xyz_ = @(n) unique(model.elems(obj_tet_(n),:));

fascicle_list = model.object_name(strncmp(model.object_name,'Fasc',4));


if ~exist('v_extracellular','var')
  %%
  
  v_extracellular = zeros(size(model.nodes(:,1))); 
  the_colormap = [G(9); (magma(99) + 0.1)/1.1];
  
  for f = fascicle_list
    this = eval(f{1});
    v_extracellular(this.idx,1) = this.pot(elec_id,:);
    
    if 0, clf %#ok<UNRCH> % debug 
      scatter3(model.nodes(this.idx,1), model.nodes(this.idx,2), ...
               model.nodes(this.idx,3), [], this.pot(elec_id,:));
    end
    
  end  
  elec_id = 1; 
  set_cs = @(c) [0.03 1]*max(abs(c));

else the_colormap = interp1([-1;0;1],[C(1,:);G(8.5);C(2,:)],linspace(-1,1,100)');
  set_cs = @(c) [-1 1]*max(abs(c));
end

if multi_geometry_model, v_extracellular = v_extracellular{1}; end

%% Generate transverse section 

z_index = 'auto'; % 1.875;
if any(named('-z')), z_index = get_('-z'); end
if ischar(z_index), 

    if strncmpi(z_index,'pos',3)
       v_max = max(v_extracellular(obj_xyz_('Fascicle1'),elec_id)); 
    elseif strncmpi(z_index,'neg',3)
         v_max = min(v_extracellular(obj_xyz_('Fascicle1'),elec_id));
    else v_max = max(abs(v_extracellular(obj_xyz_('Fascicle1'),elec_id))); 
    end

  [~,idx] = min(abs(v_max - v_extracellular(:,elec_id))); 
  z_index = model.nodes(idx,1) + 0;
end

xy_ = @(d,i) reshape(model.nodes(i,d),[],size(i,2));
cut = ~( all(xy_(1,model.elems)>z_index,2) | ...
         all(xy_(1,model.elems)<z_index,2) );
tri = [model.elems(cut,[1 2 4]); model.elems(cut,[1 2 3]); ...
       model.elems(cut,[1 3 4]); model.elems(cut,[4 3 2])]; 

cut = sum( xy_(1,tri) <= z_index, 2) == 0;
tri = double(tri(cut,:));


roi = 'auto'; % 0.2; 
if any(named('-roi')), roi = get_('-roi'); end
if ischar(roi),        roi = [0 nan nan];   
  for f = fascicle_list 
    
    x = xy_(3,obj_xyz_(f{1}));
    y = xy_(2,obj_xyz_(f{1}));
    roi = [max([roi(1); abs(x)]) nanmin([roi(2); y]) nanmax([roi(3); y])];    
  end
  
  roi = roi * [1.2 0 0; 0 1.2 -0.2; 0 -0.2 1.2];
end

switch numel(roi)
    case 1, roi = [-1 1 -0.5 1.5]*roi;
    case 2, roi = [-1 1 -0.5 1.5].*roi([1 1 2 2]);
    case 3, roi = [-1 1 1 1].*roi([1 1:3]);
    case 4, roi = roi; %#ok<ASGSL>
    otherwise roi = [-1 1 -0.5 1.5].*roi([1 1 end end]);
end

if any(named('-lin')),       si_mode = {'linear','nearest'};
elseif any(named('-patch')), si_mode = {'nearest','nearest'};
else                         si_mode = {'natural','nearest'}; 
end

avg_vxc = scatteredInterpolant(xy_(3,tri(:)),xy_(2,tri(:)), ...
                               v_extracellular(tri(:),elec_id),si_mode{:});  
[gx,gy] = meshgrid(linspace(roi(1),roi(2),500), ...
                   linspace(roi(3),roi(4),500));

clf, set(gca,'Position',[5 5 21 90]/100), hold on
imagesc(gx(1,:),gy(:,1),avg_vxc(gx,gy)), axis image xy, axis(axis)

triplot(tri, model.nodes(:,3), model.nodes(:,2),'Color',[0 0 0 0.3])
colormap(gca, the_colormap)
caxis(set_cs(caxis))
% end
%%

if any(named('-full')) 
    %%
    cut = ~( all(xy_(3,model.elems)>0,2) | all(xy_(3,model.elems)<0,2) );

    tet = model.elems(cut,:);

    tri = [model.elems(cut,[1 2 4]); model.elems(cut,[1 2 3]); ...
           model.elems(cut,[1 3 4]); model.elems(cut,[4 3 2])]; 
    cut = sum( xy_(3,tri) <= 0, 2) == 0;
    tri = double(tri(cut,:));

    % cut = ~( all(xy_(1,tri(:,1:2)) > z_max,2) | all(xy_(1,tri(:,1:2)) < z_max,2) );
    % nvc = (xy_(3,tri(:,1))-xy_(3,tri(:,2))).*(xy_(2,tri(:,1))-xy_(2,tri(:,3))) - ...
    %       (xy_(2,tri(:,1))-xy_(2,tri(:,2))).*(xy_(3,tri(:,1))-xy_(1,tri(:,3)));
    % tri = tri(nvc>0,:);
    % in_f = all(ismember(tri,obj_xyz_('Fasc')) | ismember(tri,obj_xyz_('P_F')),2);

    % trimesh(tri, m.nodes(:,3),m.nodes(:,2),m.nodes(:,1),v_meas(:,1),'FaceColor','none')
    % axis equal off, axis(axis/5)

    avg_vxc = scatteredInterpolant(xy_(1,tet(:)),xy_(2,tet(:)), ... 
                              v_extracellular(tet(:),elec_id),si_mode{:});
    [gx,gy] = meshgrid(linspace(-3,3,1000), linspace(-1,1,800));

    subplot(1,4,[2 4]), cla, hold on
    imagesc(gx(1,:),gy(:,1),avg_vxc(gx,gy)), axis image xy, axis(axis)

    triplot(tri, model.nodes(:,1),model.nodes(:,2),'Color',[0 0 0 0.3])
    colormap(gca, the_colormap) 
    caxis(set_cs(0.65))
    
elseif any(named('-mesh-3d'))
    %%
    cla, hold on
    trisurf(tri, model.nodes(:,3),model.nodes(:,2), ... 
                  v_extracellular(:,elec_id), ...
                   'EdgeColor',[0 0 0],'EdgeAlpha',0.3,'FaceAlpha',0.2)
    axis equal, axis([-0.2 0.2 -0.11 + [0 0.4]])

    for f = model.object_name(strncmp(model.object_name,'Fasc',4))

      in_f = sum(ismember(tri,obj_xyz_(f{1})),2) > 1;
      trisurf(tri(in_f,:), model.nodes(:,3),model.nodes(:,2), ... 
              model.nodes(:,1),v_extracellular(:,elec_id),'EdgeColor',[0 0 0])
    end

    in_f = all(ismember(tri,obj_xyz_('Fasc')) | ismember(tri,obj_xyz_('P_F')),2);
    lim = quantile(v_extracellular(tri(in_f,:),1),[0 1]);
    % 
    tools.tidyPlot

    caxis([0 1.1] * max(abs(lim)))
    colormap(gca, interp1([0;1],[G(8.5);C(2,:)],linspace(0,1,100)'))
    
else
  %%
  y_offset = 0;
  if any(named('-y')), y_offset = get_('-z'); 
     error TODO_add_y_offset_mark
  end
  
  x_range = 3;
  if any(named('-x')), x_range = get_('-x'); end

  
  subplot(1,4,[2 4]), cla reset

  dY = 0.15; y0 = 0; 

  for f = fascicle_list 

    x = xy_(1,obj_xyz_(f{1}));
    y = xy_(2,obj_xyz_(f{1}));

    % yc = median(y);
    yc = (max(y)+min(y))/2; 

    y = y - yc + y_offset;  y0 = y0 + dY;  
    avg_vxc = scatteredInterpolant(x,y,v_extracellular(obj_xyz_(f{1}),elec_id),si_mode{:});
    [gx,gy] = meshgrid(linspace(-x_range,x_range,401), linspace(min(y),max(y),101));

    p = avg_vxc(gx,gy); 
    % p = p - median(p(:));
    imagesc(gx(1,:),gy(:,1)+y0+yc-dY,p), hold on

    % cla, hold on
    % scatter(x,y,[],v_meas(obj_xyz_(f{1}),1),'o','filled')

    % in_f = sum(ismember(tri,obj_xyz_(f{1})),2) > 1;
    % trisurf(tri(in_f,:), m.nodes(:,3),m.nodes(:,2),m.nodes(:,1),v_meas(:,1),'EdgeColor',[0 0 0])
  end

  axis([-x_range x_range -dY y0+dY/2+yc/2]), axis xy
  tools.tidyPlot

  % ok = (all(reshape(model.nodes(model.elems,3),[],4) >= 0,2)); 
  colormap(gca, the_colormap)
  caxis(set_cs(caxis))
  % set(gcf,'Renderer','painters')

  set(gca,'YColor','none')
  plot([1 1]*min(xlim),[0 0.1],'Color',G(2),'LineWidth',2,'Clipping','off')
  text(min(xlim)-0.025,0.05,'100 µm','Color',G(2),'Rotation',90,'Horiz','center','vert','bottom')
  
  plot([1 1]*z_index, [-0.5*dY y0],':','Color',G(3),'LineWidth',1.2)

  for ee = 1:nE

    p = [min(xy_(1,model.electrode(ee).nodes)) -0.13 ...
         range(xy_(1,model.electrode(ee).nodes)) 0.03];
    rectangle('Position',p,'EdgeColor','none','FaceColor',G(3))
  end

  [~,nom] = fileparts(filename); 
  title(strrep(nom,'_','\_'),'fontweight','normal')
  
end

ch = colorbar; 
ch.Position(1) = 0.94;
ch.Parent.Children(end).CLim = caxis;



function make_this_PDF(varargin)

named = @(v) strncmpi(v,varargin,length(v)); 
get_ = @(v) varargin{find(named(v))+1};
varargin(named('-pdf')) = []; 

do_PDF  = ~any(named('-preview'));
plots.PDF_tools('setup',do_PDF); 

%%
f_ = @(x) [x.folder filesep x.name];

if any(named('-pdf-list')), list = get_('-pdf-list'); 
elseif any(named('-sens')) 
     list = dir(tools.file('sub~/eidors/sensitivity*.mat'));
else list = dir(tools.file('sub~/eidors/stimulus*.mat'));
end

%%

nRows = 5; 
if any(named('-rows')), nRows = get_('-rows'); end
  
f = figure(1); clf
set(gcf,'Position',[50 -20 760 925]) % ??? 

for ff = 1:numel(list)
  
  row = nRows - mod(ff-1,nRows);
  
  if row == nRows, f = figure(1); clf, end  
  
  figure(2); clf
  set(gcf,'Position',[820 685 760 230])
  plots.fascicle_potential_section(f_(list(ff)),varargin{:})

  h = copyobj(get(gcf,'Children'),f); 

  for ii = 1:size(h)

    p = h(ii).Position([2 4]) ./ [1 nRows]; 
    h(ii).Position([2 4]) = p.*0.75 + [row-0.8 0]*p(2);
  end


  if do_PDF && (row == 1 || ff == numel(list))
    
    close(2)    
    plots.PDF_tools('page',gcf,'page-%03d.ps',ff);    
    pause(0.1), close(gcf),  printInfo(); 
    
  else pause(0.1)
  end
end

plots.PDF_tools('compile',do_PDF,'Sensitivity_comparison (%d).pdf')




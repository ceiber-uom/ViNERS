

function s = preview_fascicles (varargin)
% preview display for mesh.insert_gmsh_fascicles
% TODO: add visualisation for 3d fascicles

named = @(v) strncmpi(v,varargin,length(v)); 

s = mesh.insert_gmsh_fascicles(varargin{:},'-info');

if isfield(s,'splines'), s = s.splines; end

if numel(s) > 1,
  s_extended = s; 
  s = s(1); 
end

%%
cla reset, hold on

nF = size(s.coeffs,3);
C = lines(max(nF,7)); G = @(v) [v v v]/10;
do_warn = []; 

if exist('s_extended','var')
  for jj = numel(s_extended):-1:2
    se = s_extended(jj);
    for ii = 1:size(se.coeffs,3)
      fill(se.outline(:,1,ii), se.outline(:,2,ii),1-G(jj),'edgeColor','none','FaceAlpha',0.5)
    end
  end
end

for ii = 1:nF
    plot(s.coeffs(1,:,ii), s.coeffs(2,:,ii),'o','Color',C(ii,:),'MarkerSize',4)
    plot(s.outline(:,1,ii), s.outline(:,2,ii),'-','Color',C(ii,:))    
    
    if any(s.outline(:,2,ii) < 0), do_warn = [do_warn ii]; end %#ok<AGROW>
end

axis equal, tools.tidyPlot, grid on, ax = axis; 



if any(named('-no-w')) || any(named('-now')), return, end

fill(ax([1 2 2 1 1]),[0 0 ax([3 3]) 0],G(3),'FaceAlpha',0.3,'EdgeColor',G(3),'LineWidth',1.1)

if ~isempty(do_warn)
    do_warn = sprintf(',#%d',do_warn);
    text(mean(xlim),mean(ylim),sprintf('Negative Y in fascicles: %s', ...
                               do_warn(2:end)),'Color','r','FontSize',14)
end

if nargout == 0, clear, end
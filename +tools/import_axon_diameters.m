
function [axons, xml_file, axon_file] = make_axons_from_EM_tracing(filename, varargin)
% [axons, xml_file, axon_file] = make_axons_from_EM_tracing(
%                                         filename_myelinated, 
%                                        [filename_unmyelinated], [...])
% 
% This code assumes length units in the input file are µm. 
% 
% Options: 
%  -opts [options] : pass options to models.axon_population for generation
%  -xml [filename] : use specified XML file for fascicle outline (disables
%                       XML generation from axon diameters) 
%  -out [filename] : set output axon file name (default: overwrite the
%                       most recent models.axon_population output)
%                    for XML generation: 
%  -e-min [0.2000] : minimum offset from hull of axons to facicle boundary
%  -e-typ [1.8634] : average offset from boundary (in µm)
%  -e-typ [0.8337] : std. deviation of offset from boundary
%                    for AXON assignment: 
%  -ratio [A C]    : what fraction of A and C fibres are afferents? 
%                     (default: 0.74 for cervical vagus nerve)
%                    aff/eff are assigned randomly unless the input table
%                    has a column "is_afferent" 
% 
% v0.1 - 2 July 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(v) strncmpi(v,varargin,length(v));
get_ = @(v) varargin{find(named(v))+1};

if nargin == 0 % automatically look in ../data for the tables

    list = dir('../data/*.xlsx'); 
    p_ = @(x) [x(1).folder filesep x(1).name]; % path expander
    is_unmy = contains(lower({list.name}),'unmy'); 
    A_fn = p_(list(find(~is_unmy,1)));
    C_fn = p_(list(find( is_unmy,1)));

elseif nargin > 1 && exist(varargin{1},'file') % 
     A_fn = filename; C_fn = varargin{1};
elseif any(named('-unmy')), A_fn = filename; C_fn = get_('-unmy');
elseif contains(filename,'_Unmy')
     A_fn = strrep(filename,'_Unmy','_My'); C_fn = filename;
else C_fn = strrep(filename,'_My','_Unmy'); A_fn = filename;
end

fprintf('Reading %s [myelinated]\n',   A_fn),  A = readtable(A_fn);
fprintf('Reading %s [unmyelinated]\n', C_fn),  C = readtable(C_fn);

A.is_axon = strcmp(A.Name,'Axon');
A.axon_id = 0*A.Area;
ok = find(~A.is_axon); 

for ii = find(A.is_axon)'

    [~,idx] = min( (A.CentroidX(ii) - A.CentroidX(ok)).^2 + ...
                   (A.CentroidY(ii) - A.CentroidY(ok)).^2 );
    A.axon_id([ii ok(idx)]) = ii; 
    ok(idx) = [];
end

all_xy = [A.CentroidX(~A.is_axon) A.CentroidY(~A.is_axon); ...
          C.CentroidX C.CentroidY];
all_od = [A.MinFeret(~A.is_axon); C.MinFeret];

nerve_outline = boundary(all_xy(:,1), all_xy(:,2));

min_excess  = 0.2    ; % All length units here µm
typ_excess = 1.8634 ; % from a random pelvic nerve fascicle
std_excess  = 0.8337 ; % (see below)

if any(named('-e-min')), min_excess = get_('-e-min'); end
if any(named('-e-typ')), typ_excess = get_('-e-typ'); end
if any(named('-e-std')), std_excess = get_('-e-std'); end

fascicle_xy = all_xy(nerve_outline,:);

u_ = @(x) x./sqrt(sum(x.^2,2));

for ii = 1:size(fascicle_xy,1)-1

    if ii == 1, this_xy = fascicle_xy([end-1 1 2],:);
    else        this_xy = [this_xy(2,:); fascicle_xy(ii+(0:1),:)];
    end % because fascicle_xy is updated inside loop

    % plot(this_xy(:,1), this_xy(:,2),'k-')
    % axis image

    vec = u_(this_xy([1 3],:) - this_xy(2,:));
    cp = cross([vec(1,:) 0],[vec(2,:) 0]);
    ov = u_( mean(vec) * sign(cp(3)));

    od = max(randn(1) * std_excess + typ_excess, min_excess); 
    od = od + all_od(nerve_outline(ii));

    fascicle_xy(ii,:) = fascicle_xy(ii,:) + od * ov;

    % hold on, plot(this_xy(2,1) + [0 ov(1)], ...
    %               this_xy(2,2) + [0 ov(2)])
   
end

fascicle_xy(end,:) = fascicle_xy(1,:);
fascicle_xy(any(isnan(fascicle_xy),2),:) = []; 

%%

smooth = @(x) conv2(x([1 1:end end],:), 1./[4;2;4],'valid');

fascicle_xy = smooth(fascicle_xy);
fascicle_xy([1 end],:) = [1;1] * mean(fascicle_xy([1 end],:)); 

clf, hold on
plot(all_xy(nerve_outline,1), all_xy(nerve_outline,2))
plot(fascicle_xy(:,1), fascicle_xy(:,2))

%% 

if any(named('-xml')), xml_file = get_('-xml');
else xml_file = strrep(strrep(A_fn,'_Myelinated_','_'),'.xlsx','.xml');
end

if ~any(named('-xml'))
    % Next steps(1) : write this XML file
    % Next steps(2) : update an axons.mat file
    fprintf('Generating %s\n', xml_file)
    f = fopen(xml_file,'wt');
    c = onCleanup(@() fclose(f));
    
    % Ripped directly from the example vagus fascicle
    
    fprintf(f, '<?xml version="1.0" encoding="ISO-8859-1"?>\n');
    fprintf(f, '<mbf version="4.0" xmlns="http://www.mbfbioscience.com" xmlns:nl="http://www.mbfbioscience.com" appname="Neurolucida 360" appversion="2020.1.1 (64-bit)">\n');
    fprintf(f, '<description><![CDATA[]]></description>\n');
    fprintf(f, '<contour name="Inner perineurium of vagus nerve" color="#FF0000" closed="true" shape="Contour">');
    fprintf(f, '  <point x="%0.2f" y="%0.2f" z="0.00" d="0.73"/>\n', fascicle_xy(2:end,:)');
    fprintf(f, '</contour>\n');
    
    % This is just placeholder, I fixed the minor bug requiring this in my own
    % version of the code but hopefully you can skip that in your code. 
    
    fprintf(f, '<contour name="Another Structure" color="#CCCCCC" closed="false" shape="Contour">\n');
    fprintf(f, '  <point x="%0.2f" y="%0.2f" z="0.00" d="0.73"/>\n', 0:3);
    fprintf(f, '</contour>\n');
    fprintf(f, '</mbf>\n');
end
%% Generate an axons file the ususal way before overwriting it with reality

% as a side effect this will always generate into a new file because 

if any(named('-op')), opts = get_('-op'); 
else opts = struct; 
end

opts.nerve.file = xml_file;
if ~isfield(opts,'axons'), opts.axons.source = 'rat-vagus-cervical'; end

mesh.insert_gmsh_fascicles('-setup',opts);
plots.preview_fascicles('-no-warning'); % Visualise
pause(0.01)

models.axon_population(opts)
axon_file = tools.file('get','axons~\ax*.mat','newest'); 
ax = load(axon_file); 

%% Figure out how positions may or may not have been transformed

close all

clf, hold on, axis equal
plot(fascicle_xy(:,1)/1e3, fascicle_xy(:, 2)/1e3)
plot(ax.nerve.outline(:,1), ax.nerve.outline(:,2) )

% the "default" processing done by mesh.insert_gmsh_fascicles is to flip
% the Y axis from the MBF default (negative, IDK why) and, if applicable,
% scale the nerve down. Then we may or may not have requested additional
% translation/rotation. Get that here as an aff

const_row = ones(size(ax.nerve.coeffs(1,:))); 

if ~isfield(ax.nerve,'source'), affine_xform = [1 0 0; 0 1 0; 0 0 1];
else affine_xform = [ax.nerve.source/1e3; const_row]' \ ...
                    [ax.nerve.coeffs; const_row]';
end

fxy_transformed = fascicle_xy/1e3.*[1 -1];
fxy_transformed(:,3) = 1; 
fxy_transformed = fxy_transformed * affine_xform; 

plot( fxy_transformed(:,1), fxy_transformed(:,2), '-k')
tools.tidyPlot, box on

% the black line should sit exactly on top of the red line. 

%% Overwrite axons.population data 

n_axons = cellfun(@numel, {ax.pop.fibre_diam}); 

AE_myel_ratio = n_axons(1)/(n_axons(1)+n_axons(2));
AE_unmy_ratio = n_axons(3)/(n_axons(3)+n_axons(4));

if any(named('-ratio')), AE_myel_ratio = get_('-ratio');
  if numel(AE_myel_ratio) > 1, 
       AE_unmy_ratio = AE_myel_ratio(2);
       AE_myel_ratio = AE_myel_ratio(1);
  else AE_unmy_ratio = AE_myel_ratio(1);
  end
end

if ~any(strcmp(A.Properties.VariableNames,'is_afferent'))
    A.is_afferent = (rand(size(A.CentroidX)) <= AE_myel_ratio); 
end
if ~any(strcmp(C.Properties.VariableNames,'is_afferent'))
    C.is_afferent = (rand(size(C.CentroidX)) <= AE_unmy_ratio);
end

for pp = 1:4, this = ax.pop(pp); % for each axon population

    % safety check
    assert(all(this.fascicle == 1),'TODO: implement multi fascicle')
    nG = max(this.size_sample); % # groups

    if this.myelinated % A fibre
        
        axon_xy = [A.CentroidX A.CentroidY 0*A.CentroidY+1]; 
        sel = A.is_axon & (A.is_afferent == this.afferent); 
        
        % if axons and fibres are not 1:1 this will likely break ... 
        ax_ids = A.axon_id(sel);         
        fi_ids = arrayfun(@(a) find( A.axon_id./~A.is_axon == a ,1 ), ax_ids); 

        this.fibre_diam = A.MinFeret(fi_ids); 
        this.axon_diam = A.MinFeret( sel );
        this.g_ratio = this.axon_diam ./ this.fibre_diam;

        % recompute sampling for axon models based on recommended approach
        gdata = [this.fibre_diam this.axon_diam]; 
        kmeans_opts = statset('MaxIter',1000);
        gm_fun = @(x,k)kmeans(x, k, 'replicate',5,'Options',kmeans_opts);    
        this.size_sample = gm_fun(gdata,nG);        
    
    else % C fibre

        axon_xy = [C.CentroidX C.CentroidY 0*C.CentroidY+1]; 
        sel = C.is_afferent == this.afferent;

        this.fibre_diam = C.MinFeret(sel);
        diam_qtiles = quantile(C.MinFeret, 2*nG+1); 
        diam_qtiles(1:2:end) = [];

        [~,this.size_sample] = arrayfun(@(d) min(abs(diam_qtiles - d)), ... 
                                            this.fibre_diam);
    end

    axon_xy = axon_xy(sel,:) .* [1e-3 -1e-3 1]; 
    axon_xy = axon_xy * affine_xform; 

    this.axon_xy = axon_xy(:,1:2);
    this.fascicle = ones(size(this.fibre_diam));

    plot( axon_xy(:,1), axon_xy(:,2), '.','Color', this.color)

    ax.pop(pp) = this; 
end

ax.opts = [ax.opts {'Axon population data imported from:', ...
                    A_fn, C_fn}]; 

if any(named('-out')), axon_file = get_('-out'); end
if exist(axon_file,'file'), verb = 'Overwrite'; else verb = 'Save'; end
fprintf('%sing %s ...\n', verb(1:end-1), axon_file); 
save(axon_file,'-struct','ax'); 

if nargout == 0, clear
else axons = ax;
end

return



function script_get_outline_offset_prior %#ok<DEFNU> 
%% Load a random Pelvic Nerve fascicle to get a prior for the offset

clear
R19_401_folder = 'E:/Microscopy Data/R19-401-L fascicle2';
p_ = @(x) [x(1).folder filesep x(1).name]; % path expander

file_ = @(s) p_(dir([R19_401_folder '/'  s]));

A = readtable(file_('*_myelinated_*.xlsx'));
C = readtable(file_('*_unmyelinated_*.xlsx'));

all_xy = [A.CentroidX A.CentroidY; C.CentroidX C.CentroidY];

all_od = [A.FiberDiameter_1; C.MinFeretDistance];
nerve_outline = boundary(all_xy(:,1), all_xy(:,2));
nerve_outline(1) = []; 

%%
filename = file_('*.xml');

if ~exist('xml','var')
    if isempty(which('tools.file')), path(path, 'E:/ViNERS/code'), end
    fprintf('Loading %s\n ... ', filename)
    xml = tools.parse_xml(filename); 
    fprintf('Done!\n')
end

% lambda functions to extract values
the = @(x,n) x.Children(strncmpi({x.Children.Name},n,length(n)));
attr_ = @(x,n) x.Attributes(strcmpi(n,{x.Attributes.Name})).Value;

contours = the(xml,'contour');
c_names = arrayfun(@(c) attr_(c,'name'), contours, 'UniformOutput', false);
is_fascicle = strcmp( c_names, 'Fascicle area' );
perinuerium = contours(is_fascicle); 
perinuerium = the(perinuerium,'point');

peri_xy = [str2double(arrayfun(@(p) attr_(p,'x'), perinuerium,'unif',0));
           str2double(arrayfun(@(p) attr_(p,'y'), perinuerium,'unif',0))]';

clear the attr_ contours c_names is_fascicle perinuerium

%%

offset_dist = @(p) sqrt(min(sum((all_xy(p,:) - peri_xy).^2,2)));
offset_um = arrayfun(offset_dist, nerve_outline);

excess_offset = offset_um - all_od(nerve_outline);

plot(offset_um); 
hold on, plot(all_od(nerve_outline))

%%
sel = (1:60); sel(excess_offset(sel)<0) = []; 

clf
histogram( excess_offset(sel), 12 )

[~,p] = vartest2( excess_offset(sel), offset_um(sel) );

mean(excess_offset(sel))
std(excess_offset(sel))

%%
ex = 1; 
est = fitdist(excess_offset(sel) .^ ex , 'normal');
qqplot( excess_offset(sel) .^ ex , est)







% this example shows how to set up and run a ViNERS experiment

clear
tools.file('set','sub-1') % what subject? 
% tools.file('set','sub-1','sam-2') % example subject + sample

%%

% read example configuration file
opts = tools.parse_json(tools.file('~/source/example.json')); 

% You could also set up your model configuration with something like this:
% 
% setup = struct; 
% setup.array = models.electrode_array;
% setup.mesh.MeshLengthMax = 0.1 % maximum mesh size, mm
% setup.mesh.DomainSize = [6 6 3 0.3] % x y z and thickness of tissue
% setup.nerve.file = '~/source/nerve/rat-pelvic-nerve.xml
% ... 

if strcmp(tools.configuration('name'),'THE-STREAM')    
    opts.all.no_parallel = true; % issue on my home PC. 
end

%% Send this configuration to mesh.insert_gmsh_fascicles. 

mesh.insert_gmsh_fascicles('-setup',opts);
plots.preview_fascicles % Visualise

opts.nerve.Perineurium_um = mesh.get_perineurium_thickness; % read from file

%% Compute Ve / Tm 
% Validated 13 June 2021 CDE [commit 23c8dee]

if ~tools.model_status('nerve_anatomy','done')
    tools.cache('reset') % reset cache of temporary files
    models.nerve_anatomy(opts) % run models.nerve_anatomy to generate mesh and fields
    % the output of this step is saved to tools.file('eidors~\')
    
    %% Visualise sensitivity function 
    clf, plots.view_sensitivity('-3d');
    caxis([0 max(caxis)])
end

%% Generate axon population for model
% Validated 13 June 2021 CDE [commit e74229d]

if ~tools.model_status('axon_population','done')
    models.axon_population(opts)
end

% the output of this step is saved to tools.file('axons~\')
% 
% by default, the axon population is named based on the nerve of origin 
% for the axon sample (e.g. pelvic-nerve). This can be changed e.g. : 
% models.axon_population(opts,'-out','alternate-sample')

% A nicer visualisation of the model axon population can be generated
% using:
% 
% >> plots.preview_axons

%% Generate representitive membrane current profiles for axons
% Validated 15 June 2021 CDE [commit f9660c4] -no-parallel

if ~tools.model_status('membrane_currents','done')
    models.membrane_currents(opts);
end

%% Determine stimulation thresholds for each axon. 
% as above, test fully when parallel toolbox debugged.  

if ~tools.model_status('axon_thresholds','done')
    stim_file = tools.file('get','eidors~\stim*.mat'); 
    models.axon_thresholds(opts, '-file',stim_file)
end

% The output of this step is saved to tools.file('thresholds~\')
% 
% models.axon_thresholds uses a default stimulus of a 100 µs pulsewidth
% square biphasic bipolar stimulus on E12 (however that is defined on your
% array) with a 50 µs interphase gap. For each axon, this stimulus is
% modulated up and down (logarithmic search) until a 1% increase in
% stimulus pushes the axon superthreshold; the greater of these two values
% is reported as the threshold

%% Determine stimulation responses for each axon to fixed stimuli.
% as above, test fully when parallel toolbox debugged.  

if ~tools.model_status('nerve_stimulation','done')    
    stim_file = tools.file('get','eidors~\stim*.mat'); 
    models.nerve_stimulation(opts,'-file',stim_file)
end

% The output of this step is saved to tools.file('stim~\').
% 
% models.nerve_stimulation differs from models.axon thresholds in that it
% uses fixed, pre-determined levels of stimulation as would be done in a
% typical in-vivo experiment. Thresholds are estimated but are not as exact
% as those estimated from models.axon_thresholds. 

%% Detemine single-fibre action potential properties for each axon.
% Validated 15 June 2021 CDE [commit 5900258] -no-parallel

if ~tools.model_status('axon_sfap','done')    
    sens_file = tools.file('get','eidors~\sens*.mat'); 
    models.axon_sfap(opts, '-file', sens_file)    
end

% the output of this step is saved to tools.file('sfap~\')

%% Simulate ensemble recordings of the activity in many axons.
% Validated 15 June 2021 CDE [commit 7b6c2e2] -no-parallel

if ~tools.model_status('nerve_recording','flat','1','done')
    sens_file = tools.file('get','eidors~\sens*.mat');     
    models.nerve_recording(opts, '-file', sens_file)
end
 
% the output of this step is saved to tools.file('waves~\')

% the waves saved by models.nerve_recording are seperated by axon class 
% and fascicle of origin; this allows post-hoc composition of recordings
% under different combinations of activity. 

% models.compose_waves assembles the output of models.nerve_recording 
%   (or models.ECAP_response) into the final form as needed by plots.

% plots.preview_waves suggusts that this isn't working correctly yet even
% if it technically runs without crashes. [ ] TODO continue this 


%% Simulate recordings of electrically-evoked compound action potentials.

if ~tools.model_status('ECAP_recording','done')
    sens_file = tools.file('get','eidors~\sens*.mat');
    models.ECAP_recording(opts, '-file', sens_file)
end

% the output of this step is also saved to tools.file('waves~\')


%%

% visualise your results 





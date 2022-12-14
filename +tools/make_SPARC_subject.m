
function make_SPARC_subject(number, varargin)

if nargin == 0, root = tools.file('sub~/','sub=new'); 
elseif ischar(number)
  root = sprintf('%s%s/',tools.file('~/primary/'),number); 
  if nargin > 1
    sample_run = sprintf('%s/',varargin{:});
    root = [root sample_run];
  end
else root = sprintf('%s-%d/',tools.file('~/primary/sub'),number);    
end

if exist([root 'thresholds/'],'dir'), return, end

disp('Generating structure for:')
disp(root)

mkdir([root 'axons/'])
mkdir([root 'eidors/'])
mkdir([root 'waves/'])
mkdir([root 'stim/'])
mkdir([root 'thresholds/'])

w_ = @(varargin) fprintf(varargin{:});

f = fopen([root 'axons/README.md'],'wt'); 
w_(f,'This folder contains the following:\n');
w_(f,'- axons.mat     : output of make_axon_population for my particular');
w_(f,                 ' axon distribution and fascicles\n');
w_(f,'- /MRG/*.mat    : output by models.membrane_currents\n');
w_(f,'- /Sundt/*.mat  : output by models.membrane_currents\n');
w_(f,'- /Gaines/*.mat : output by models.membrane_currents\n');
w_(f,'\nInformation current as of 05-June-2020');
w_(f,'\nv0.1  Calvin Eiber  ceiber@unimelb.edu.au');
fclose(f);

f = fopen([root 'eidors/README.md'],'wt'); 
w_(f,'This folder contains the following:\n');
w_(f,'- sensitivity (n).mat  : output of models.nerve_anatomy\n');
w_(f,'- stimulus (n).mat     : output of models.nerve_anatomy -stimulus\n');
w_(f,'- elec-geom (n).mat    : optional, descripton of array geometry');
w_(f,                        ' input to models.nerve_anatomy');
w_(f,'\nInformation current as of 05-June-2020');
w_(f,'\nv0.1  Calvin Eiber  ceiber@unimelb.edu.au');
fclose(f);

f = fopen([root 'waves/README.md'],'wt'); 
w_(f,'This folder contains simulated population responses, grouped by');
w_(f,' simulation type, generated by models.ECAP_response or models.ongoing_activity.\n\n');
w_(f,'One or more of the following may be present:\n');
w_(f,' - /base/epoch_b00_c00 (n).mat : baseline (typ. 0.1 Hz)\n');
w_(f,' - /sweep/epoch_k00_c00 (n).mat : firing rate sweep (typ. 0.1-20 Hz)\n');
w_(f,' - /default/epoch_k00_c00 (n).mat : default firing rates (typ. 0.1,1,10 Hz)\n');
w_(f,' - /burst/epoch_k00_f00_c00 (n).mat : non-stationary pop. firing rate, bursts of different size/frequency\n');
w_(f,' - /step/epoch_b00_k00_c00_w00 (n).mat : non-stationary pop. firing rate, fast and slow step changes\n');
w_(f,' - /stim/stim_CL000.mat : different simulated ECAP responses. If "CL" notation was used to define current levels, this syntax is used');
w_(f,' - /stim/stim_UA000.mat : different simulated ECAP responses. "UA" is the stimulus amplitude in microamps\n');

w_(f,'\nepoch filenames contain some clues as to their contents:\n');
w_(f,'- _b: baseline firing rate\n');
w_(f,'- _k: peak firing rate (or average stationary firing rate)\n');
w_(f,'- _c: within-population coherence\n');
w_(f,'- _f: bursting frequency\n');
w_(f,'- _w: non-stationary response slope\n');

w_(f,'\nThese definitions can be changed by user-controlled software so look to the documentation for individual experiments if there is something that disagrees with the information presented here.\n')

w_(f,'\nInformation current as of 15-November-2021');
w_(f,'\nv0.2  Calvin Eiber  ceiber@unimelb.edu.au');
fclose(f);




if 1
  %%
  f = fopen([root 'thresholds/README.md'],'wt'); %%#ok<UNRCH>
  w_(f,'This folder contains the following (output by models.axon_thresholds)\n');
  w_(f,'in folders named \\stimulus (stimulus-metadata)\\ following the contents of eidors~\\\n:');
  w_(f,'- MRG-fascicle{{1-n}}.mat    : output \n');
  w_(f,'- Sundt-fascicle{{1-n}}.mat  : output \n');
  w_(f,'- Gaines-fascicle{{1-n}}.mat : output \n');
  w_(f,'\nInformation current as of 15-Jan-2021');
  w_(f,'\nv0.2  Calvin Eiber  ceiber@unimelb.edu.au\n');
  
  w_(f,'\nNOTE: this is different to the results generated by models.nerve_stimulation, which\n'); 
  w_(f,'  still estimate threshold but use an approach which mirrors that of an in-vivo experiment\n');
  w_(f,'models.membrane_currents -stimulus is an old synonym for models.axon_thresholds.\n');
  
  fclose(f);
end

if 1
  %%
  f = fopen([root 'stim/README.md'],'wt');   
  w_(f,'This folder contains the following (output by models.nerve_stimulation)\n');
  w_(f,'in folders named \\stimulus (stimulus-metadata)\\ following the contents of eidors~\\\n:');
  w_(f,'- MRG-fascicle{{1-n}}.mat    : output. Each files contains spiketime responses to all stimuli. \n');
  w_(f,'- Sundt-fascicle{{1-n}}.mat  : output \n');
  w_(f,'- Gaines-fascicle{{1-n}}.mat : output \n');
  w_(f,'\nNOTE: this is different to the results generated by models.axon_thresholds, which \n'); 
  w_(f,'  compute each axon''s individual threshold. These files mirror the structure of an in-vivo\n');
  w_(f,'  experiment, using a series of pre-set current levels.');
  w_(f,'\nNOTE: as of this readme (rep) is not embedded in the filename here or in the downstream\n');
  w_(f,'  waves~/ folder, using different stimulus configurations will require different folder names.\n');
  w_(f,'\nInformation current as of 15-Jan-2021');
  w_(f,'\nv0.2  Calvin Eiber  ceiber@unimelb.edu.au\n');
  
  fclose(f);
end

% axon_model
% membrane_currents 

f = fopen([root 'INPUTS'],'wt'); 
w_(f,'# This file captures the model inputs for this run\n');
w_(f,'# This file is used by tools.INPUT_file(list,key)\n');
w_(f,'\n# v0 Generated on %s',datestr(now));
fclose(f);

% [x] models.axon_model < called by models.membrane_current
% [ ] models.membrane_current    < outputs /axons, /thresholds
% [ ] models.pelvic_nerve        < outputs /eidors
% [ ] models.population_response < outputs /waves
% [x] models.random_raster < called by models.population_response
% [x] models.random_raster < called by models.population_response


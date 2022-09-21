
% This script fixes the EIDORS overloads issue if Andy Adler and the eidors
% gang haven't gotten around to fixing this themselves yet. 
% 
% This script modifies the EIDORS source code! use with caution ... once. 
% 
% Step 1: all of the .m files in the eidors/overloads folder get "eidors_"
% prepended to the filename. 
% 
% Step 2: all of the eidors code gets searched for references to the munged
% overloads. any line of code which uses a modified file is changed to
% point to the new file.
% 
% Step 3: profit.

s_ = @(varargin) system(sprintf(varargin{:})); % system call
p_ = @(x,varargin) [x.folder filesep x.name varargin{:}]; % path expander

cd(tools.configuration('eidors')); 

list_to_change = dir('./overloads/*.m'); 

for ii = 1:numel(list_to_change)            
    new_name{ii} = ['eidors_' list_to_change(ii).name];   
    if contains(new_name{ii},'eidors_eidors_'), % avoid double changes
        new_name{ii} = ''; continue, 
    end
    s_('REN "%s" "%s"', p_(list_to_change(ii)), new_name{ii});    
end

new_name = strrep(new_name,'.m','\('); 
old_name = strrep(new_name,'eidors_','');

%% This section 

list_to_search = dir('./**/*.m'); 

for ii = 1:numel(list_to_search)
    
    mfile = fopen(p_(list_to_search(ii)),'r+t'); 
    code = fread(mfile,inf,'*char')'; 
    frewind(mfile); 
    
    ncc = numel(code);    
    for jj = 1:numel(new_name)
        code = regexprep(code,['([^\w])' old_name{jj}],['$1' new_name{jj}]);
    end
    
    if numel(code) == ncc, fclose(mfile); continue, end
    
    n_hits = (numel(code)-ncc) / 7 ; % numel('eidors_')
    fprintf('%s: made %d replacements\n', list_to_search(ii).name, n_hits)
    fprintf(mfile,'%s',code);
    fclose(mfile);
end
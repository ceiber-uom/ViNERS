
function value = isOctave

  persistent cached_value % memoise
  if ~isempty(cached_value), value = cached_value; return, end

  value = str2num(tools.configuration('octave')); %#ok<ST2NM>
  if isempty(value), value = false; end % Default: NOT octave
  cached_value = value;
  
  return
  
  %% This is the other way of finding this nugget of information out

  if isdeployed   %#ok<UNRCH>
       value = isempty(strfind(ctfroot,'MATLAB')) %#ok<STREMP>
  else value = isempty(which('MATLAB'));   
  end
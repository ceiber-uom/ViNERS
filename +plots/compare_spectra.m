

function compare_spectra(varargin)
% plots.compare_spectra(...)
% appears to be a duplicate of plots.compare_recordings
% (as in plots.compare_recording absorbed this functionality)
% 
% CDE 18 Nov 2021

if nargin == 0, varargin = {'-burst','-quick'}; end

named = @(v) strncmpi(v,varargin,length(v)); 

if evalin('caller','exist(''pop_info'',''var'')')
  pop_info = evalin('caller','pop_info');
  pop_response = evalin('caller','pop_response');
else
  plots.compare_recordings(varargin{:},'-spectra','-gather-data');
  assignin('caller','pop_info',pop_info)         %#ok<NODEF>
  assignin('caller','pop_response',pop_response) %#ok<NODEF>
end
%%

do_pdf = any(named('-pdf')); 
do_normalised = any(named('-norm'));
plots.PDF_tools('setup',do_pdf); 

for i_type = 1:size(pop_response,1)  % For each set of active axons 

  hz = pop_response(1).hz;
  
  nC = pop_info.page.nRows; 
  nR = max(pop_info.page.id); 

  figure, clf
  for rr  = 1:nR % row/column structure defined by plots.response_waves
    
    this_row = find(pop_info.page.id == rr)';     
    
    for cc = 1:nC
      %% Make each row on the image
      idx = this_row(pop_info.page.get_row_(cc));
      colors = pop_info.page.get_color_(cc); 

      subplot(nR,nC,cc + (rr-1)*nC); cla, hold on
      
      for ii = idx
      
        sel = (pop_info.group == ii); 
        
        spectra = [pop_response(i_type,sel).spect]; 
        
        if ii == idx(1), yRef = mean(spectra,2); end
        
        if do_normalised          
          spectra = spectra ./ yRef;
        end
        
        
        sd = std(spectra,[],2) ; % /sqrt(size(spectra,2)); 
        y0 = mean(spectra,2);        
        sd(isnan(sd)) = nanmedian(sd(:));
        
        ys_neg =  flipud(y0-sd); 
        ys_neg(ys_neg < 0) = min(y0); 
        
        fill([hz fliplr(hz)], [y0+sd;ys_neg],colors(idx == ii,:), ... 
                      'EdgeColor','none','FaceAlpha',0.3)
        plot(hz,y0,'-','Color',colors(idx == ii,:),'LineWidth',1.5)
        
        
      end
      
      set(gca,'YScale','log','userdata','primary')      
      axis tight, tools.tidyPlot
      
      p = get(gca,'Position')./[1 1 3 4];
      axes('Position',p+[2 0.1 0 0].*p([3 4 3 4]));
      k = gausswin(7); k = k./sum(k);
      hold on
      for ii = idx
        
        sel = (pop_info.group == ii); 

        stime = pop_response(1).spk_time;
        srate = cat(1,pop_response(i_type,sel).spk_rate); 
        srate = mean(srate,1);
        srate = conv([srate fliplr(srate) srate],k,'same'); 
        srate = fliplr(srate(numel(stime) + (1:numel(stime))));
        
        
        plot(stime,srate,'-','Color',colors(idx == ii,:),'LineWidth',1.1);
      end
      set(gca,'userdata','inset')
      axis tight, tools.tidyPlot
      set(gca,'XTick','','YTick','')
      
    end    
  end


  h = findobj(gcf,'UserData','primary'); 
  set(h,'YLim',[min([h.YLim]) max([h.YLim])])
  set(h,'YTick',10.^(-10:10))

  sel = sum(cat(1,h.Position),2) == min(sum(cat(1,h.Position),2));
  set(h(~sel),'XTickLabel','','YTickLAbel','')    
  xlabel(h(sel),'Frequency (hz)'), 
  if do_normalised, ylabel(h(sel),'fold change'), lbl = 'Normalised Spectra';
  else                ylabel(h(sel),'power (a.u.)'), lbl = 'Power Spectra';
  end    
  h = findobj(gcf,'UserData','inset'); 
  set(h,'YLim',[min([h.YLim]) max([h.YLim])])

  % sel = sum(cat(1,h.Position),2) == min(sum(cat(1,h.Position),2));
  
  txt = sprintf(', %s', pop_info.axon_type{pop_info.active_ty{i_type}});  
  txt = sprintf('%s: %s', txt(3:end), lbl);
  suptitle(txt)


  if do_pdf
       plots.PDF_tools(gcf,'a%03d-summary.ps',i_type);
       pause(0.05), close(gcf)
  else pause(0.05)
  end
end % active_ID

if do_pdf
  plots.PDF_tools('combine','get','Population %s (%%d).pdf','next', lbl)
end

%%
pro drawing_2d,Um,pres,vx,vz, time, itt , dt

@nonspec_c.pro
@grid_c
@init_c

if firsttime then begin
    ;if it's the first time we draw we make sure everything is well configured
    
    if !d.name eq 'WIN' then begin

        ; ...................................................
        ; (a) making sure the system is only using 256 colors
        device,decomposed=0
        ; ...................................................

        ; ...................................................
        ; (b) loading a certain palette, i.e., color table (standard 
        ;       color tables are, e.g.,  39, or 1 or 2) 
        loadct,39,/silent
        ;     (to check out all available color tables, 
        ;               say in the idl prompt: xloadct)
        ; ...................................................
        ; window size conserving the aspect ratio of the fisical domain:
        winxsiz= 1200d0 * (xf-x0)/(zf-z0)   &    winzsiz= 1200d0
        ;compute the plot magnification to coincide with the window size
        pltmag = [(winxsiz-40)/(2*npx) , (winzsiz-40)/(2*npz)]
        ; ...................................................
        ; (c) opening windows if not available yet:
        device,window_state=wstat
        ; opening the windows that we are going to use (1 in this case)
        for iwin=0,0 do if wstat[iwin] eq 0 then window,iwin,xsiz=winxsiz,ysiz=winzsiz,title='plots'
        ; ...................................................
        ; ...................................................
        ; (d)change global settings for the plot variables
        !p.charsize = 2 & !p.charthick=1.5 & & !x.charsize=1.5 & !y.charsize=1.5
        !p.thick = 3 & !x.thick=3 & !y.thick=3
        !P.multi    = [0, 2, 2]
        !p.background = 255
        !p.color  = 0
        !x.margin = [8,1] & !y.margin = [5,2]
        ; ...................................................

        ; When leaving this if-block, set 'firsttime' to zero so that 
        ;      it is never entered again during the execution of the code. 
        firsttime = 0
    endif else begin
        print, 'drawing.pro: there is a problem in the !d.name variable'
        print, 'unable to draw properly'
        firsttime = 0
    endelse
endif

;computing the curves to draw normalized perturbed ones
umnorm = (Um-Um00)/Um00 
pnorm  = (pres-p00)/p00
vxnorm  = vx/cs00
vznorm  = vz/cs00

;the amplituds must be equal except for a gamma factor 
;so we compute one initial amplitud to use it as a frame
ampini = max(abs((Uminit - Um00)/Um00))
;the initial velocity makes the y axes in the velocity to move
;towards another equilibrium position [v0-amp,v0+amp] so we compute them
v_ampini = vxinit/cs00
med   = (max(v_ampini) + min(v_ampini))/2
range =  max(v_ampini) - min(v_ampini)
    
;we create the frames and plot the densities, velocities and pressures
if plttype eq 'contour' then begin
  
  contour, umnorm,xx,zz,/isotropic ,title='(!4 q-q!3!l00!n)/!4q!3!l00'     ,ytitle='z',xtickformat='(A1)',$
    nlev=200,c_colors = [50:250:1],/yst,/xst, xr = [x0,xf] ,yr = [z0,zf],/fill
  if itype eq 'packet' and track eq 'yes' then oplot, [xc],[zc], psym=6, SYMSIZE=2
  if itype eq 'packet'                    then oplot, xtrac[0:-20],ztrac[0:-20]
  ;------------------------------------------------------------------------------------------------------
  contour, pnorm ,xx,zz,/isotropic ,title='(p-p!l00!n)/p!l00!n' , ytickformat='(A1)'  ,xtickformat='(A1)',$
    nlev=200,c_colors = [50:250:1],/yst,/xst, xr = [x0,xf] ,yr = [z0,zf],/fill
  if itype eq 'packet' and track eq 'yes' then oplot, [xc],[zc], psym=6, SYMSIZE=2
  if itype eq 'packet'                    then oplot, xtrac[0:-20],ztrac[0:-20]
  ;------------------------------------------------------------------------------------------------------
  contour, vxnorm,xx,zz,/isotropic ,title='v!lx!n/C!lS00!n'     ,ytitle='z',xtitle='x',$
    nlev=200,c_colors = [50:250:1],/yst,/xst, xr = [x0,xf] ,yr = [z0,zf],/fill
  if itype eq 'packet' and track eq 'yes' then oplot, [xc],[zc], psym=6, SYMSIZE=2
  if itype eq 'packet'                    then oplot, xtrac[0:-20],ztrac[0:-20]
  ;------------------------------------------------------------------------------------------------------
  contour, vznorm,xx,zz,/isotropic ,title='v!lz!n/C!lS00!n'     ,xtitle='x'           , ytickformat='(A1)',$
    nlev=200,c_colors = [50:250:1],/yst,/xst, xr = [x0,xf] ,yr = [z0,zf],/fill
  if itype eq 'packet' and track eq 'yes' then oplot, [xc],[zc], psym=6, SYMSIZE=2
  if itype eq 'packet'                    then oplot, xtrac[0:-20],ztrac[0:-20]
  ;------------------------------------------------------------------------------------------------------
 
  xyouts,/norm,0.45,0.53 ,CHARSIZE = 3, 'itt= '  + strtrim(string(itt),2)
  xyouts,/norm,0.45,0.49, CHARSIZE = 3, 'time= ' + strtrim(string(time,form='(f7.2)'),2)
endif else if plttype eq 'color' then begin
  
  tvscl, congrid(umnorm,npx*pltmag[0],npz*pltmag[1]) , 15               ,30+npz*pltmag[1]
  tvscl, congrid(pnorm ,npx*pltmag[0],npz*pltmag[1]) , 30+npx*pltmag[0] ,30+npz*pltmag[1]
  tvscl, congrid(vxnorm,npx*pltmag[0],npz*pltmag[1]) , 15               ,15
  tvscl, congrid(vznorm,npx*pltmag[0],npz*pltmag[1]) , 30+npx*pltmag[0] ,15
  
  xyouts,/norm,0.53,0.53 , CHARSIZE = 2, 'itt= '  + strtrim(string(itt),2)
  xyouts,/norm,0.53,0.47, CHARSIZE = 2, 'time= ' + strtrim(string(time,form='(f7.2)'),2)
endif else begin
  print, 'drawing.pro : "' + plttype + '" Type of plot not suported, NOT PLOTTING ANITHING'
endelse

wait,0.005

end




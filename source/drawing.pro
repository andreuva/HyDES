pro drawing,Um,pres,vv, grid, time,itt

@nonspec_c.pro
@grid_c
@init_c

if firsttime then begin
    ;if it's the first time we draw we make sure everything is well configured
    set_plot, 'WIN'
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

        ; ...................................................
        ; (c) opening windows if not available yet:
        device,window_state=wstat
        ; opening the windows that we are going to use (1 in this case)
        for iwin=0,0 do if wstat[iwin] eq 0 then window,iwin,xsiz=1200,ysiz=1200,title='plots'
        ; ...................................................
         
        ; ...................................................
        ; (d)change global settings for the plot variables
        !p.charsize= 3 & !x.charsize=1.5 & !y.charsize=1.5
        !P.multi = [0, 1, 3]
        !p.background = 255
        !p.color=0
        !p.thick = 3 & !x.thick=3 & !y.thick=3  & !p.charthick=1.5
        !x.margin = [18,5] & !y.margin = [8,2]
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
vnorm  = vv/cs00

;the amplituds must be equal except for a gamma factor 
;so we compute one initial amplitudE to use it as a frame
if plttype eq 'cutz' then ampinim = max(abs((Uminit[npx/2,*]  - Um00)/Um00))
if plttype eq 'cutx' then ampinim = max(abs((Uminit[*,npz/2]  - Um00)/Um00))
ampinip = max(abs((presinit - p00)/p00))
;the initial velocity makes the y axes in the velocity  move
;towards another equilibrium position [v0-amp,v0+amp] so we compute them
vvinit = sqrt(vxinit*vxinit + vzinit*vzinit)
v_ampini = vvinit/cs00
med   = (max(v_ampini) + min(v_ampini))/2
range = max(v_ampini) - min(v_ampini)

axes = 'Z'
if plttype eq 'cutx' then begin
  axes = 'X'
  zc = -1d18
endif
if track ne 'yes' then zc = -1d18

;we create the frames and plot the densities, velocities and pressures
;--------------------------------------------------------------------------------
plot ,grid,vnorm       ,/nodata ,yr=[med - range*0.80 , med + range*1.6],/yst,$
  xr = [min(grid),max(grid)],/xst,  title= 'VELOCITY', ytit = 'v/C!lS00!n',$
  xtickformat='(A1)',ymar=[2,2]; , xtit = 'Z'
oplot,grid,vnorm       ,color=85
oplot,[min(grid),max(grid)],[med - range*0.5 , med - range*0.5],line=3, thick=0.7
oplot,[zc,zc],[-1d8, 1d8]                             ,line=3, thick=0.7
oplot,[min(grid),max(grid)],[med + range*0.5 , med + range*0.5],line=3, thick=0.7
;--------------------------------------------------------------------------------
plot ,grid,Umnorm       ,/nodata, yr=[-2*ampinim, +2*ampinim],/yst, $
   xr = [min(grid),max(grid)],/xst, title= 'DENSITY', ytit = '(!4 q-q!3!l00!n)/!4q!3!l00',$;  xtit = 'Z'
  xtickformat='(A1)',ymar=[2,2]
oplot,grid,Umnorm       ,color=150
oplot,[min(grid),max(grid)],[+ampinim,+ampinim] ,line=3,thick=0.7
oplot,[min(grid),max(grid)],[0, 0]            ,line=3,thick=0.7
oplot,[zc,zc],[-1d8, 1d8]                     ,line=3, thick=0.7  
oplot,[min(grid),max(grid)],[-ampinim,-ampinim] ,line=3,thick=0.7
;--------------------------------------------------------------------------------
plot ,grid,pnorm        ,/nodata, yr=[-1.3*ampinip, +1.3*ampinip],/yst,$
   xr = [min(grid),max(grid)],/xst, title= 'PRESSURE', ytit = '(p-p!l00!n)/p!l00!n',xtit = axes 
oplot,grid,pnorm        ,color=190
oplot,[min(grid),max(grid)],[+ampinip,+ampinip] ,line=3,thick=0.7
oplot,[min(grid),max(grid)],[0, 0]                      ,line=3,thick=0.7
oplot,[zc,zc],[-1d8, 1d8]                               ,line=3, thick=0.7
oplot,[min(grid),max(grid)],[-ampinip,-ampinip] ,line=3,thick=0.7
;--------------------------------------------------------------------------------
; write the ITERATION and the time
xyouts,/norm,0.8,0.65 , CHARSIZE = 3, 'itt= '  + strtrim(string(itt),2)
xyouts,/norm,0.78,0.32, CHARSIZE = 3, 'time= ' + strtrim(string(time,form='(f7.2)'),2)
wait,0.005

end




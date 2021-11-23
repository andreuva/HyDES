; pro MAIN
;------------------------------------------------------------------
; 2D HYDRODINAMICAL SIMULATOR   || AUTHOR: ANDRES VICENTE AREVALO
;------------------------------------------------------------------
; -------------
; COMMON BLOCKS 
; -------------
@nonspec_c
@grid_c
@init_c
;-----------------
; COMPILE RUTINES
;-----------------
@drawing
@drawing_2d
@initrout_2d
@advance_fg_2d
@deriv2d
@boundc_2d
@midval
@fluxes
@primitives
@tstep
@rain
@ray_tracing
; ------------------
; READ IN PARAMETERS
; ------------------
; OPEN THE INPUT DATA FILE
filename = 'parameter_files/datain_2d.dat'
openr,unit,/get_lun,filename
strd=''   ; declare strd as a string variable for later use

; READ IN THE PARAMETERS FOR THE NUMERICAL GRID and CREATE IT
readf,unit,nintx  & npx=nintx+2d0        
readf,unit,nintz  & npz=nintz+2d0
readf,unit,x0,xf
readf,unit,z0,zf  

; READ IN MAX ITERATIONS, MAX TIME, OUTPUT PERIODICITY AND PLOT TYPE
readf, unit, itmax     
readf, unit, timef 
readf, unit, plot_cad, store_cad
readf, unit,form='(a15)',strd & plttype   = strtrim(strd,2)

readf, unit, strd ; just jump over an empty line in the input data file

; READ IN THE PARAMETERS OF THE INITIAL AND BOUNDARY CONDITION 
readf, unit,form='(a15)',strd & itype     = strtrim(strd,2)
readf, unit,form='(a15)',strd & shape     = strtrim(strd,2)
readf, unit,form='(a15)',strd & track     = strtrim(strd,2)
readf, unit,form='(a15)',strd & boundcond = strtrim(strd,2)
readf, unit,ondx ,ondz   
readf, unit,wp ,wt 
readf, unit,xc ,zc
readf, unit,nwave   
readf, unit,ampl    
readf, unit,v0_cs

;parameters of the simulated 'rain' flag and cadence
readf, unit,form='(a15)',strd & rain_flag = strtrim(strd,2)
readf, unit,rain_cad

readf, unit, strd ; jump over an empty line

; READ IN THE COURANT-FRIEDRICHS-LEWY PARAMETER
readf, unit, cflparam

free_lun,unit ; close and deallocate the unit
; --------------------------------------
; READ PARAMETERS FROM INPUT FILE -- END
; --------------------------------------
; --------------------------------------
; CALCULATE THE GRID
; --------------------------------------
; COMPUTE THE GRID AND THE INCREMENT IN ZZ
dx = (xf-x0)/nintx
dz = (zf-z0)/nintz

xx = dindgen(npx)*dx + (x0 - (dx/2d0))
zz = dindgen(npz)*dz + (z0 - (dz/2d0))

gridx = xx#(dblarr(npz)+1)
gridz = (dblarr(npx)+1)#zz
; --------------------------------------
; CALCULATE THE GRID -- END
; --------------------------------------
; ------------------------------------------------
; CHECK IF THE INITIAL PARAMETERS ARE PROPERLY SET
; ------------------------------------------------
if (zf le z0) or (xf le x0) then begin
  print, 'Error in datain.dat: Initial points are larger than lasts ones'
  print, 'Z0 = ', z0, ' zf = ', zf
  print, 'x0 = ', x0, ' xf = ', xf
  stop
endif

if (nintx le 1) or (nintz le 1) then begin
  print, 'Error in datain.dat: The Number of intervals should be 2 or greater'
  print, 'Nintx = ',nintx
  print, 'Nintz = ',nintz
  stop
endif

; make sure that the wavenumbers are integers and are lees than the lenght
ondx = ulong64(ondx)
ondz = ulong64(ondz)

if (ondx gt nintx) or (ondz gt nintz) then begin
  print, 'Error in datain.dat: The Number of wave N (ond(x,z))) in k = 2pi*ond(x,z)*/L should be Nint(x,z) or smaller'
  print, 'ondx = ',ondx,'  nintx = ',nintx   &   print, 'Nintz = ',ondz,'  nintz = ',nintz
  print, 'ondx and ondz has been redefined as nintx and nintz.'
  ondx = nintx      &      ondz = nintz 
endif

;check the boundary conditions and puts the 'periodic' as default
if (boundcond ne 'periodic') and (boundcond ne '0deriv') then begin
  print, 'Boundari conditions "'+ boundcond +'" not valid, computing perodic ones as default.'
  boundcond = 'periodic'
endif

;check for a warning if the stability can be compromise
if cflparam ge 0.99d0 then print, 'WARNING: Cfl param is "',cflparam,'" stability may be compromised'

; warn if no wave is selected in the wave modes
if (ondx eq 0) and (ondz eq 0) and (itype eq 'sound wave' or itype eq 'chladni' or itype eq 'packet') then begin
  print, 'WARNING: There are no modes selected kx=0,kz=0 : kz=1 has been set'
  ondz=1
endif

;warn if the rain is activated outside a gaussian mode or a incompatible plot mode is selected
if itype eq 'packet' then boundcond = '0deriv'
if rain_flag eq 'yes' and itype ne 'gaussian' then print, 'WARNING: RAIN IS ACTIVATED IN A NON GAUSSIAN MODE'
; ------------------------------------------------------
; CHECK IF THE INITIAL PARAMETERS ARE PROPERLY SET - END
; ------------------------------------------------------
; -------------------------
;  INITIAL CONDITIONS 
; -------------------------
firsttime=1

initrout   ;the routine initrout sets up the initial condition.

; this calculates the initial values of the 'densities' for momentum and energy:
Upxinit  = Uminit*vxinit
Upzinit  = Uminit*vzinit 
Ueinit = presinit/(gamm-1d0) + Uminit*(vxinit^2d0+vzinit^2d0)/2d0 

; this calculates the sound speed array for the initial condition:
csinit = sqrt(gamm*presinit/Uminit)

;draw and store the initial conditions
if plttype eq 'cutz' then begin
  vv = sqrt(vxinit*vxinit + vzinit*vzinit)
  drawing,Uminit[fix(npx/2),*],presinit[fix(npx/2),*],vv[fix(npx/2),*], zz, 0.0, 0.0
endif else if plttype eq 'cutx' then begin
  vv = sqrt(vxinit*vxinit + vzinit*vzinit)
  drawing,Uminit[*,fix(npz/2)],presinit[*,fix(npz/2)],vv[*,fix(npz/2)], xx, 0.0, 0.0
endif else begin
  drawing_2d,Uminit,presinit,vxinit,vzinit,0.0,0.0,0.0
endelse
save,/VARIABLES, FILENAME = 'results/states/inital_conditions.sav'
WRITE_PNG, 'results/plots/initial_conditions.png', TVRD(/TRUE)
; -------------------------
;  INITIAL CONDITIONS - end
; -------------------------
; ------------------------------------------------------------
; DECLARE and/or INITIALIZE VARIABLES JUST BEFORE STARTING THE LOOP
; ------------------------------------------------------------
vx    = vxinit  & vz   = vzinit & cs   = csinit   & pres = presinit  
Um    = Uminit  & Ue   = Ueinit & Upx  = Upxinit  & Upz  = Upzinit
Umn   = Um      & Uen  = Ue     & Upxn = Upx      & Upzn = Upz 

; the 'densities', in the sense of the conservation laws, 
;      are 'Um', 'Upx', 'Upz' and  'Ue'. 

itt  = 0L   
time = 0d0 
; --------------------------------------------------------
; INITIALIZE VARIABLES JUST BEFORE STARTING THE LOOP - END
; --------------------------------------------------------
; ---------------
; BIG LOOP BEGINS
; ---------------
while itt lt itmax and time lt timef do begin   
    
;   CALCULATE FLUXES FROM DENSITIES ('fluxes' is the name of the subroutine) 
    fluxes,Um,Upx,Upz,Ue,    fmx,fmz,   fpxx,fpxz,fpzz,fpzx,  fex,fez
    
;   TIMESTEP 
    dt=tstep(Um,Upx,Upz,cs)

    if (time+dt) gt timef then begin
      dt = timef-time
    endif

;   ADVANCE ONE TIMESTEP using the chosen numerical scheme 
    advance, dt, Um,Upx,Upz,Ue, fmx,fmz,   fpxx,fpxz,fpzz,fpzx,  fex,fez, Umn,Upxn,Upzn,Uen

;   BOUNDARY CONDITIONS 
    boundc,boundcond,'left',Umn    & boundc,boundcond,'right',Umn
    boundc,boundcond,'left',Uen    & boundc,boundcond,'right',Uen
    boundc,boundcond,'up',Umn      & boundc,boundcond,'down',Umn
    boundc,boundcond,'up',Uen      & boundc,boundcond,'down',Uen
    
    boundc,boundcond,'left',Upxn   & boundc,boundcond,'right',Upxn
    boundc,boundcond,'left',Upzn   & boundc,boundcond,'right',Upzn
    boundc,boundcond,'up',Upxn     & boundc,boundcond,'down',Upxn
    boundc,boundcond,'up',Upzn     & boundc,boundcond,'down',Upzn       

;   CALCULATE PRIMITIVE VARIABLES FROM THE DENSITIES    
    primitives, Umn,Upxn,Upzn,Uen, vxn,vzn,presn

;   EXCHANGE NEW AND OLD VARIABLES
    Um = Umn & Upx = Upxn & Upz = Upzn & Ue = Uen & pres = presn  & vx = vxn  & vz = vzn 
    
    itt++      ; advance itterations
    time += dt ; advance time
    
    ;evolve the central point of the packet to track it
    if itype eq 'packet' then begin
      indz = (zc-z0)/dz
      zc = zc + interpolate(cs,npx/2,indz)*dt
    endif
    
;   PLOT RESULTS 
    if itt mod plot_cad eq 0 or itt eq 0 or time ge timef then begin       
        if plttype eq 'cutz' then begin
          vv = sqrt(vx*vx + vz*vz)
          drawing,Um[fix(npx/2),*],pres[fix(npx/2),*],vv[fix(npx/2),*], zz, time, itt
        endif else if plttype eq 'cutx' then begin
          vv = sqrt(vx*vx + vz*vz)
          drawing,Um[*,fix(npz/2)],pres[*,fix(npz/2)],vv[*,fix(npz/2)], xx, time, itt
        endif else begin
          drawing_2d,Um,pres,vx,vz, time, itt, dt
        endelse
    endif
    
;   STORE RESULTS - to store results in a file from time to time
    if itt mod store_cad eq 0 then begin
      save,/VARIABLES, FILENAME = 'results/states/estate_time_' $
        + strtrim(string(time,form='(f7.2)'),2)+'.sav'
      WRITE_PNG, 'results/plots/plot_time'+strtrim(string(itt+1d4,form='(f7.1)'),2)+'.png', TVRD(/TRUE)
    endif
    
    ;make rain if the mode is selected:)
    rain, rain_flag, rain_cad, itt, Um,pres,vx,vz   ,Upx,Upz,Ue,cs
    
endwhile
; --------------
; BIG LOOP - end
; --------------
print,'program finished.'+$
      '    itt= '+ strtrim(string(itt),2) +$
      '; time = '+ strtrim(string(time),2)
end

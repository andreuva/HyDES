; ----------------------------------------------------------------------
; ROUTINE EXACT_SOLUTION
;
;    PURPOSE: Calculate the exact solution of a sound wave.
;
;    INPUT ARGUMENTS: time
;    
;    COMMON BLOCKS: some inputs are taken from common blocks as:
;    xf,z0,zz from grid 
;    cs00,ampl, itype and shape from init_c
;    gamma from nonspec
;
;    OUTPUT:  the fluxes fmz, fpz, fez
; ----------------------------------------------------------------------

pro exact_solution, time, ex,ex_v,zline,line

@nonspec_c.pro
@grid_c
@init_c
  
if itype eq 'sound wave' then begin
  ;compute the analitical solution for a cosine sound wave
  kk = 2d0*!dpi/(zf-z0)
  phy= 2d0*!pi/5d0 -  kk*z0 + 2*kk*zz

  if shape eq 'cosine' then begin
    ;in the first mode the velocity is in phase with the other cuantities
    ex   = cos(-(v00+cs00)*kk*time - kk*zz + phy)
    ex_v = v00/cs00 + ampl*cos(-(v00+cs00)*kk*time - kk*zz + phy)

    ;computing the line in the maximum moving with the sound speed
    zline = (v00+cs00)*time*[1,1] + (z0 - 2d0*!dpi/(5d0*kk)) + 4d0*!dpi/kk
    ;apling boundary conditions
    zline = (zline-z0)mod(zf-z0) + z0
    ;setting the hight of the line
    line = [-1d3,1d3]

  endif else if shape eq 'second mode' then begin
    ;in the second mode the velociti is out of phase by pi
    ex  = cos(-(v00-cs00)*kk*time - kk*zz + phy)
    ex_v= v00/cs00 + ampl*cos((-v00+cs00)*kk*time - kk*zz + phy + !dpi)

    ;computing the line in the maximum moving with the sound speed
    zline = (v00-cs00)*time*[1,1] + (z0 - 2*!dpi/(5d0*kk)) + 0d0*!dpi/kk
    ;apling boundary conditions
    zline = (zline - zf)mod(zf-z0) + zf
    line  = [-1d3,1d3]

  endif else begin
    ;taking into account other shapes (althought we did it in introut routine)
    ;just in case :)
    print, 'drawing.pro: Shape of sound wave not supported.'
    stop
  endelse
  ;if it's not a sound wave we can't compute the solution so we initialice the variables
  ;to 0 just for avoid errors in the plots
endif else if itype eq 'pulse' then begin
  ;computing the exact solution as a initial conditión moving by the sound speed
  
  ;calculate the amount travel by the pulse and aplying boundary conditions
  zp = zz - time*(v00 + cs00) 
  zp = (zp-zf)mod(zf-z0) + zf
  
  ;compute the function as displace initial condition
  z1  = (5.5d0*z0 + 1.5d0*zf)/7d0
  z2  = (3d0*z0   + 4d0  *zf)/7d0
  W = 0.3d0
  hz  = Beselj((-20d0*zp/5d0),0)*( 1+tanh((zp-z1)/W) )*( 1-tanh((zp-z2)/W) )
  
  ;compute the exact solution like the function times the amplitud  
  ex    = hz
  ex_v  = v00/cs00 + ampl*hz
  
  ;seting maximum values to 0 to not plot them
  zline = [0]
  line  = [0]
  
endif else begin
  
  zline = [0]
  ex    = [0]
  ex_v  = [0]
  line  = [0]
  
endelse

end
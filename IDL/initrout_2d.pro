; ----------------------------------------------------------------------
; ROUTINE INITROUT
;
;    PURPOSE:  Calculate the arrays of density, velocity and pressure
;                 at time t=0.
;    INPUT ARGUMENTS: 
;                    ALL ARGUMENTS PASSED VIA COMMON BLOCKS
;         INIT_C COMMON:
;              itype:   char variable with the choice of initial condition
;              shape:   char variable with some subsidiary choice
;              Um00 :   value of the equilibrium density
;              V00  :   value of the equilibrium velocity 
;              v0_cs:   value in wich v00 is prop to cs00
;              p00  :   value of the equilibrium pressure
;              cs00 :   value of the speed of sound in equilibrium
;              ampl :   amplitud of the pertubation
;        GRID COMMON
;              z0, zf:  inital and final values of the range
;        NONSPEC COMMON
;              gamm :   adiabatic exponent for our gass
;
;    OUTPUT:  VIA INIT_C COMMON:
;             Uminit:  array with the initial values of the density
;             vvinit:  array with the initial values of the velocity
;             presinit:array with the initial values of the pressure
; ----------------------------------------------------------------------

PRO INITROUT

;call the common blocks
@nonspec_c
@grid_c
@init_c

; initialize some of the variables in the common blocks
gamm = 5d0/3d0  &   p00 = 1.15d0   &   Um00 = 0.8d0
cs00 = sqrt(gamm*p00/Um00)

; set the equilibrium velocity
v00 = v0_cs*cs00

;create the matrix to compute the perturbation
hz = dblarr(npx,npz)
hz_vx = dblarr(npx,npz)
hz_vz = dblarr(npx,npz)

;choose the right type of initial perturbation
if itype eq 'sound wave' then begin
  ;if its a sound wave compute the initial conditions with a cosine
  ;or if it's the second mode, a cosine + a phase in the velocity
  if shape eq 'cosine' then begin
    
    ;go through all the points in the grid to compute the perturbation
    foreach xi,xx,ii do begin
      foreach zi,zz,jj do begin
        hz[ii,jj]    = cos(2d0*!dpi*ondx*xi/(xf-x0) + 2d0*!dpi*ondz*zi/(zf-z0))
        hz_vx[ii,jj] = cos(2d0*!dpi*ondx*xi/(xf-x0) + 2d0*!dpi*ondz*zi/(zf-z0))*cos(atan(double(ondz)/double(ondx)))
        hz_vz[ii,jj] = cos(2d0*!dpi*ondx*xi/(xf-x0) + 2d0*!dpi*ondz*zi/(zf-z0))*sin(atan(double(ondz)/double(ondx)))
      endforeach
    endforeach
    
    ;compute the densities and pressures with the vale of the equilibrium + the perturbation
    Uminit   = Um00  + Um00     *ampl*hz
    presinit = p00   + gamm*p00 *ampl*hz
    
    ;if the wave is in one axis make sure we dont have the other axis perturbation
    if ondz eq 0 then begin
      vxinit   = v00   + cs00     *ampl*hz_vx
      vzinit   = 0*v00   + cs00     *ampl*hz_vz*0
    endif else if ondx eq 0 then begin
      vxinit   = 0*v00   + cs00     *ampl*hz_vx*0
      vzinit   = v00   + cs00     *ampl*hz_vz
    endif else begin
      vxinit   = v00   + cs00     *ampl*hz_vx
      vzinit   = v00   + cs00     *ampl*hz_vz
    endelse
    
  endif else if shape eq 'second mode' then begin

    foreach xi,xx,ii do begin
      foreach zi,zz,jj do begin
        hz[ii,jj]    = cos(2d0*!dpi*ondx*xi/(xf-x0) + 2d0*!dpi*ondz*zi/(zf-z0))
        hz_vx[ii,jj] = cos(2d0*!dpi*ondx*xi/(xf-x0) + 2d0*!dpi*ondz*zi/(zf-z0) + !dpi)*cos(atan(double(ondz),double(ondx)))
        hz_vz[ii,jj] = cos(2d0*!dpi*ondx*xi/(xf-x0) + 2d0*!dpi*ondz*zi/(zf-z0) + !dpi)*sin(atan(double(ondz),double(ondx)))
      endforeach
    endforeach
    
    Uminit   = Um00  + Um00     *ampl*hz
    presinit = p00   + gamm*p00 *ampl*hz
    
    if ondz eq 0 then begin
      vxinit   = v00   + cs00     *ampl*hz_vx
      vzinit   = 0*v00   + cs00   *ampl*hz_vz*0
    endif else if ondx eq 0 then begin
      vxinit   = 0*v00   + cs00   *ampl*hz_vx*0
      vzinit   = v00   + cs00     *ampl*hz_vz
    endif else begin
      vxinit   = v00   + cs00     *ampl*hz_vx
      vzinit   = v00   + cs00     *ampl*hz_vz
    endelse
    
  endif else begin
    ;add an else stament to prevent errors for other shapes that are not implemented
    print,'routine initrout: the "' + shape + '" shape for sound wave is not defined yet'
    stop
  endelse

;if gaussian is choosed a gaussian function is computed as a initial perturbation
endif else if itype eq 'gaussian' then begin

  ;compute the sigma and central point in relation with the fisical domain choosed
  sigma = mean([(xf-x0),(zf-z0)])/15d0
  xp = (xf+x0)/2d0    &     zp = (zf+z0)/2d0
  
  ;go through all the point and compute the gaussian and the gradient for the velocities
  foreach xi,xx,ii do begin
    foreach zi,zz,jj do begin
      hz[ii,jj]    = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))
      hz_vx[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))*2d0*(xi-xp)/sigma
      hz_vz[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))*2d0*(zi-zp)/sigma
    endforeach
  endforeach
  
  ;compute the initial states as eq. + perturbation
  Uminit   = Um00  + Um00     *ampl*hz
  presinit = p00   + gamm*p00 *ampl*hz
  vxinit   = v00   + cs00     *ampl*hz_vx
  vzinit   = v00   + cs00     *ampl*hz_vz
  
  ;if v00 is not set to 0 warn that the gaussian will move te center over time
  if v00 ne 0 then print, 'The gaussian will be moving in the direction of k' $
    + ' because of the inital velociti v0_cs in datain.dat'
  
endif else if itype eq 'chladni' then begin
  ;this compute a Chladni mode of vibration of our experimental space
  ;in reality there are no Chladni modes in a fluid but just for fun
  
  if ondx lt 1 then ondx=1
  if ondz lt 1 then ondz=1
  
  foreach xi,xx,ii do begin
    foreach zi,zz,jj do begin
      hz[ii,jj]    = sin(2*!dpi*ondx*xi/(xf-x0))*sin(2*!dpi*ondz*zi/(zf-z0))
      hz_vx[ii,jj] = sin(2*!dpi*ondx*xi/(xf-x0))
      hz_vz[ii,jj] = sin(2*!dpi*ondz*zi/(zf-z0))
    endforeach
  endforeach
  
  Uminit   = Um00  + Um00     *ampl*hz
  presinit = p00   + gamm*p00 *ampl*hz
  vxinit   = v00   + cs00     *ampl*hz_vx
  vzinit   = v00   + cs00     *ampl*hz_vz

;if we want to simulate a shock
endif else if itype eq 'shock' then begin
  
  ;compute the jump relations of the shock and set the Mach number
  if v0_cs lt 1 then v0_cs = 2
  Mo = v0_cs
  p1 = (2d0*gamm*Mo^2d0/(gamm+1d0) - (gamm-1d0)/(gamm+1d0))*p00
  Um01= (2d0/Mo^2d0/(gamm+1d0) + (gamm-1d0)/(gamm+1d0))*Um00
  
  ;chose over the diferent types of shocks
  if shape eq 'gaussian' then begin
  ;Shock with a gaussian profile
    sigma = mean([(xf-x0),(zf-z0)])/15
    xp = (xf+x0)/2d0
    zp = (zf+z0)/2d0
    foreach xi,xx,ii do begin
      foreach zi,zz,jj do begin
        hz[ii,jj]    = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))
        hz_vx[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))*2d0*(xi-xp)/sigma
        hz_vz[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))*2d0*(zi-zp)/sigma
      endforeach
    endforeach
  
    Uminit   = Um00  + (Um00 - Um01)*hz
    presinit = p00   + (p1   - p00)*hz
    vxinit   = Mo*cs00*hz_vx/max(hz_vx)
    vzinit   = Mo*cs00*hz_vz/max(hz_vz)
  
  endif else if shape eq 'table gaus' then begin
  ;shock with a table function with gaussian wings.
    sigma = 1/sqrt(2d0*alog(2d0))
    xp = (xf+x0)/2d0
    zp = (zf+z0)/2d0
    foreach xi,xx,ii do begin
      foreach zi,zz,jj do begin
        
        if sqrt((xi-xp)*(xi-xp) + (zi-zp)*(zi-zp)) lt 1 then begin
          hz[ii,jj]    = 1
          hz_vx[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0)/2d0)*(xi-xp)/sigma
          hz_vz[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0)/2d0)*(zi-zp)/sigma
        endif else begin
          hz[ii,jj]    = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0)/2d0)*2d0
          hz_vx[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0)/2d0)*(xi-xp)/sigma
          hz_vz[ii,jj] = exp(-(((xi-xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0)/2d0)*(zi-zp)/sigma
        endelse
        
      endforeach
    endforeach

    Uminit   = Um00  + (Um00 - Um01)*hz
    presinit = p00   + (p1   - p00)*hz
    vxinit   = Mo*cs00*hz_vx/max(hz_vx)
    vzinit   = Mo*cs00*hz_vz/max(hz_vz)
  
  endif else if shape eq 'table tanh' then begin
    ;shock with a tanh() function as a table function
    xp = (xf+x0)/2d0
    zp = (zf+z0)/2d0
    foreach xi,xx,ii do begin
      foreach zi,zz,jj do begin
        rr = sqrt((xi-xp)^2d0 + (zi-zp)^2d0)
        hz[ii,jj]    =   (1d0-tanh((rr-1d0)/0.2d0 ))/2d0
        hz_vx[ii,jj] =   (1d0-tanh((rr-1d0)/0.2d0 ))/2d0*rr*cos(atan((zi-zp),(xi-xp)))
        hz_vz[ii,jj] =   (1d0-tanh((rr-1d0)/0.2d0 ))/2d0*rr*sin(atan((zi-zp),(xi-xp)))
      endforeach
    endforeach

    Uminit   = Um00  + (Um00 - Um01)*hz
    presinit = p00   + (p1   - p00)*hz
    vxinit   = Mo*cs00*hz_vx/max(abs(hz_vx))
    vzinit   = Mo*cs00*hz_vz/max(abs(hz_vz))    
    
  endif else begin
    ;add an else stament to prevent errors for other shapes that are not implemented
    print,'routine initrout: the "' + shape + '" shape for shock is not implemented yet (gaussian,table gaus or table tanh)'
    stop
  endelse

;if we want to send a wave packet
endif else if itype eq 'packet' then begin
  
  ;create the matrix to compute the discontinuities in the medium to see the packet travel
  disc  = dblarr(npx,npz)   &     disc_2= dblarr(npx,npz)
  xc    = x0 + (xf-x0)*xc    &     zc = z0 + (zf-z0)*zc
  
  ;compute the wave vector and its angle
  k0 = 2*!dpi*nwave/sqrt((xf-x0)^2d0 + (zf-z0)^2d0)
  theta = atan(double(ondz),double(ondx))
  
  ;set some initial parameters of the wave and the medium
  ;wt  = 1d9         &     wp  = 1d0
  Um00= 1.1d0       &     Um01= 0.05d0
  z1  =-2           &      z2 = 4.0d0
  AA = (1/sqrt(Um00) - 1/sqrt(Um01))/2d0   &    BB = (1/sqrt(Um00) + 1/sqrt(Um01))/2d0
  
  ;compute the perturbation and the discontinuity in the medium
  foreach xi,xx,ii do begin
    foreach zi,zz,jj do begin
      
      hz[ii,jj]    =   exp(-(((xi-xc)*cos(theta) + (zi-zc)*sin(theta)))^2d0/wp^2d0)*$
        exp(-((-(xi-xc)*sin(theta) + (zi-zc)*cos(theta)))^2d0/wt^2d0)*$
        cos(k0*((xi-xc)*cos(theta) + (zi-zc)*sin(theta)))
      hz_vx[ii,jj] =   hz[ii,jj]*cos(theta)
      hz_vz[ii,jj] =   hz[ii,jj]*sin(theta)
      
      disc[ii,jj]  = (1d0 - tanh((zi-z1)/0.2d0 ))/2d0
      
      if zi le z1 then begin
        disc_2[ii,jj]= Um00
      endif else if zi le z2 then begin
        disc_2[ii,jj]= (AA*cos(!dpi*(zi-z1)/(z2-z1)) + BB)^(-2d0)
      endif else begin
        disc_2[ii,jj]= Um01
      endelse
    endforeach
  endforeach
  
  ;chose for the tipe of discontinuity in the medium (or homogeneus by default) 
  if shape eq 'discont' then begin
    Uminit   = Um00  + Um00*ampl*hz + (-Um00 + Um01)*disc
    Um00     = Um00 + (-Um00 + Um01)*disc 
  endif else if shape eq 'gradient' then begin
    Uminit   =         Um00*ampl*hz +  disc_2
    Um00     = disc_2
  endif else begin
    Uminit   = Um00  + Um00*ampl*hz
  endelse
  
  presinit = p00   + gamm*p00 *ampl*hz
  vxinit   = v00   + cs00     *ampl*hz_vx
  vzinit   = v00   + cs00     *ampl*hz_vz
  
  ;compute the sound speed in the equilibrium and the path that will folow the wave
  csinit = sqrt(gamm*presinit/Um00)
  ray_tracing, xc , zc, k0*cos(theta), k0*sin(theta), csinit

  if plttype eq 'cutz' and (shape eq 'discont' or shape eq 'gradient') then Um00 = Um00[npx/2,*] 
  if plttype eq 'cutx' and (shape eq 'discont' or shape eq 'gradient') then Um00 = Um00[*,npz/2]

endif else begin
  ;finally we take into account other types not impleted yet to prevent errors
  print,'routine initrout: the "'+itype+'" itype is not defined yet'
  stop
endelse
end
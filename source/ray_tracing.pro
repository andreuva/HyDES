; ----------------------------------------------------------------------
; ROUTINE INITROUT
;
;    PURPOSE:  Calculate the path of a sound wave acros the medium
;    
;    INPUT ARGUMENTS:
;              
;              xinit , zinit    : Starting point of the computation
;              kxinit, kzinit   : Wave vector of the perturbation
;              csinit           : Matrix of the sound speed in the medium
;         

;        GRID COMMON
;              z0, zf:  inital and final values of the range
;              dx,dz :  values of increment in the grid
;
;    OUTPUT:  VIA PLOT FIGURE OF THE WAVE VECTORS IN THE PATH OF THE WAVE:
; ----------------------------------------------------------------------
pro ray_tracing, xinit , zinit, kxinit, kzinit, csinit

  ;call the necesary common
  @grid_c
  @nonspec_c
  
  ; initialice the variables for solving the advance equation
  xp = xinit    &     zp = zinit   
  kx = kxinit   &     kz = kzinit   &     cspdx=csinit     &     cspdz=csinit
  xs = [xp]     &     zs = [zp]     &     kxs = [kxinit]   &     kzs = [kzinit]
  
  ; compute the increment in time to use in the integration
  dt = min([dx,dz])/max(csinit)
  
  ;compute the derivates of the sound speed for doing the gradient
  deriv2d,csinit ,derx = cspdx    & deriv2d,csinit ,derz = cspdz
  
  ;iterate over the computation until the end of the domain is reached
  while xp gt x0 and xp lt xf and zp gt z0 and zp lt zf do begin
    
    ;compute the index of the point to interpolate
    indx = (xp-x0)/dx
    indz = (zp-z0)/dz
    
    ;advance the wave vector with the derivate of the sound speed
    kxn = kx - sqrt(kx*kx+kz*kz) * interpolate(cspdx,indx,indz) * dt
    kzn = kz - sqrt(kx*kx+kz*kz) * interpolate(cspdz,indx,indz) * dt
    
    ;compute the components of the unitari vector ek
    ekx = kxn/sqrt(kxn*kxn + kzn*kzn)
    ekz = kzn/sqrt(kxn*kxn + kzn*kzn)
    
    ;advance the point   
    xn  = xp + interpolate(csinit,indx,indz) * ekx * dt
    zn  = zp + interpolate(csinit,indx,indz) * ekz * dt

    ;update the wave vector and the points to the new ones
    xp = xn     &     zp = zn
    kx = kxn    &     kz = kzn
    
    ;store the result to plot it later
    kxs = [kxs,kx]   &     kzs = [kzs,kz]
    xs = [xs,xp]     &     zs = [zs,zp]
  endwhile
  ;plot the points traveled by the wave with its corresponding wave vectors
  gfx = vector(kxs[0:-1:20],kzs[0:-1:20],xs[0:-1:20],zs[0:-1:20],xrange=[x0,xf],yrange=[z0,zf], $
    title='RAY TRACING OF THE PACKET', xtitle='X', ytitle='Z' )
  ;save the trayectori to plot in the advance to see the match
  xtrac = xs      &     ztrac = zs    
end
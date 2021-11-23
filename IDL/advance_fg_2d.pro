; ----------------------------------------------------------------------
; ROUTINE ADVANCE
;
;    PURPOSE:  Calculate variables at timestep n+1.
;              THIS ROUTINE CONTAINS THE NUMERICAL SCHEME USED IN THE CODE.
;
;    INPUT ARGUMENTS:  - the densities at timestep n, namely Um, Up, Ue
;                      - the fluxes fmz, fpz, fez
;
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks (like the grid parameters and array zz, etc)
;
;    OUTPUT:  the densities at timestep n+1, namely Umn, Upn, Uen
; ----------------------------------------------------------------------

pro advance, dt, Um,Upx,Upz,Ue, fmx,fmz,   fpxx,fpxz,fpzz,fpzx,  fex,fez, Umn,Upxn,Upzn,Uen

  @nonspec_c.pro
  @grid_c

;  initialaice the arrays wich will contain the derivates of the fluxes  
  fmxdx  = fmx  & fmzdz  = fmz
  fexdx  = fex  & fezdz  = fez
  fpxxdx = fpxx & fpzzdz = fpzz 
  fpzxdx = fpzz & fpxzdz = fpxz
  
  ;calculate the derivates of the mas, energy and momentum fluxes in x and z
  deriv2d, fmx , derx = fmxdx     &     deriv2d, fmz , derz = fmzdz 
  deriv2d, fex , derx = fexdx     &     deriv2d, fez , derz = fezdz
  deriv2d, fpxx, derx = fpxxdx    &     deriv2d, fpzz, derz = fpzzdz
  deriv2d, fpzx, derx = fpzxdx    &     deriv2d, fpxz, derz = fpxzdz
  
  ;compute the densities in the time t + dt/2
   Umn[0:-2,0:-2]  = midval(Um)  - (dt/2)*(fmxdx  + fmzdz)
   Uen[0:-2,0:-2]  = midval(Ue)  - (dt/2)*(fexdx  + fezdz)
  Upxn[0:-2,0:-2]  = midval(Upx) - (dt/2)*(fpxxdx + fpxzdz)
  Upzn[0:-2,0:-2]  = midval(Upz) - (dt/2)*(fpzzdz + fpzxdx)
  
  ;calculate the fluxes in t + dt/2
  fluxes,Umn[0:-2,0:-2],Upxn[0:-2,0:-2],Upzn[0:-2,0:-2],Uen[0:-2,0:-2],  $
      fmxn,fmzn,   fpxxn,fpxzn,fpzzn,fpzxn,  fexn,fezn
  
  ;calculate the derivates of the mas, energy and momentum in x and z in t+dt/2
  deriv2d, fmxn , derx = fmxdx     &     deriv2d, fmzn , derz = fmzdz
  deriv2d, fexn , derx = fexdx     &     deriv2d, fezn , derz = fezdz
  deriv2d, fpxxn, derx = fpxxdx    &     deriv2d, fpzzn, derz = fpzzdz
  deriv2d, fpzxn, derx = fpzxdx    &     deriv2d, fpxzn, derz = fpxzdz
  
  ;compute the densities in the time t + dt 
   Umn[1:-2,1:-2] =  Um[1:-2,1:-2] - dt*(fmxdx  + fmzdz)
   Uen[1:-2,1:-2] =  Ue[1:-2,1:-2] - dt*(fexdx  + fezdz)
  Upxn[1:-2,1:-2] = Upx[1:-2,1:-2] - dt*(fpxxdx + fpxzdz)
  Upzn[1:-2,1:-2] = Upz[1:-2,1:-2] - dt*(fpzzdz + fpzxdx)

end

; ----------------------------------------------------------------------
; ROUTINE FLUXES
;
;    PURPOSE: Calculate the fluxes fmz, fpz, fez that appear in the
;             conservation laws as functions of the densities.
;
;    INPUT ARGUMENTS: the densities Um, Up and Ue
;
;    OUTPUT:  the fluxes fmz, fpz, fez
;
; ----------------------------------------------------------------------

pro fluxes,Um,Upx,Upz,Ue,    fmx,fmz,   fpxx,fpxz,fpzz,fpzx,  fex,fez
  
  ;compute the velocity and pressure for compute the fluxes
  primitives, Um,Upx,Upz,Ue,  vx,vz,pres ; this calculates vx,vz and pres

  fmx = Upx                  ; mass flux
  fmz = Upz
  
  fpxx = Upx*Upx/Um + pres    ; momentum flux
  fpxz = Upx*Upz/Um 
  fpzz = Upz*Upz/Um + pres
  fpzx = Upz*Upx/Um 
  
  fex = (Ue + pres) * vx    ; energy flux
  fez = (Ue + pres) * vz

end

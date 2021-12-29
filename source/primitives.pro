; ----------------------------------------------------------------------
; ROUTINE primitives
;
;    PURPOSE:  Calculate the primitive variables vv and pres from the
;              mass, momentum and Uey densities.
;
;    INPUT ARGUMENTS:  ARRAYS WITH THE VALUES OF:
;                     Um  : Density of matter
;                     Up  : Density of momentum
;                     Ue  : Density of energy
;
;    COMMON BLOCKS: gamm (the adiabatic exponent) is passed via common blocs
;
;    OUTPUT:
;           vv  : array with the velocities
;           pres: array with the pressures
; ----------------------------------------------------------------------
pro primitives, Um,Upx,Upz,Ue,  vx,vz,pres

@nonspec_c.pro

vx   = Upx/Um
vz   = Upz/Um

pres = (gamm-1d0)*(Ue - (Upx*Upx + Upz*Upz)/(2d0*Um))

end

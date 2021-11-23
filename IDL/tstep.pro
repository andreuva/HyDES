; ----------------------------------------------------------------------
;
; FUNCTION CFL
;
;    PURPOSE: Calculate the timestep to guarantee numerical stability
;
;    INPUT ARGUMENTS: those necessary to calculate the sound speed, namely
;           Um, Up, pres
;
;    COMMON BLOCKS: Note that some necessary input is passed via common
;                       blocks (like gamm or the grid zz)
;    OUTPUT:  the delta t
; ----------------------------------------------------------------------

FUNCTION tstep,Um,Upx,Upz,cs

@nonspec_c.pro
@grid_c.pro

;compute the characteristic velocities
vcharac1 = abs(sqrt(Upx*Upx + Upz*Upz)/Um + cs)
vcharac2 = abs(sqrt(Upx*Upx + Upz*Upz)/Um - cs)

;compute the proper delta t with the cflparam
dt = cflparam * min([dx,dz]) / max([vcharac1, vcharac2])

return, dt

end

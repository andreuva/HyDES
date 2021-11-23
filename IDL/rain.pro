; ----------------------------------------------------------------------
; ROUTINE RAIN
;
;    PURPOSE:  introduce random points of rain in our simulated fluid
;    
;    INPUT ARGUMENTS:
;         RAIN_FLAG: a flag for computing the rain or not
;         RAIN_CAD : a mesure of how often in itteratios the rain occurs
;         
;         The state of the simulation as:
;         itt, Um,pres,vx,vz   ,Upx,Upz,Ue,cs
;         
;                    ARGUMENTS PASSED VIA COMMON BLOCKS
;        INIT COMMOM
;              Um00, p00, cs00 initial values to compute the perturbation
;        GRID COMMON
;              x0, xf, z0, zf:  inital and final values of the range
;              npx, npz      :  points of the grid
;        NONSPEC COMMON
;              gamm :   adiabatic exponent for our gass
;
;    OUTPUT:  VIA INPUTS Um,pres,vx,vz   ,Upx,Upz,Ue,cs
; ----------------------------------------------------------------------

pro rain, rain_flag, rain_cad, itt, Um,pres,vx,vz   ,Upx,Upz,Ue,cs

  ;call the common blocks
  @nonspec_c
  @grid_c
  @init_c
  
  ; decide if we generate a 'drop' or not with the flag,cadence and a random choice
  if (rain_flag eq 'yes') and (itt mod rain_cad eq 0) and (randomu(!null,1) le 0.5) then begin

    
    ;chose the size and location of the new 'drop'
    tam = RANDOMU(!NULL,1)*20d0 + 10d0
    sigma = mean([(xf-x0),(zf-z0)])/tam
    xp = RANDOMU(!NULL,1)*(xf-x0)*0.9d0 + x0
    zp = RANDOMU(!NULL,1)*(zf-z0)*0.9d0 + z0
    
    ;compute the gaussian perturbation over the actual state
    hz = dblarr(npx,npz)
    hz_vx = dblarr(npx,npz)
    hz_vz = dblarr(npx,npz)
    foreach xi,xx,ii do begin
      foreach zi,zz,jj do begin
        hz[ii,jj]    = exp(-(((xi - xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))
        hz_vx[ii,jj] = exp(-(((xi - xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))*2d0*(xi - xp)/sigma
        hz_vz[ii,jj] = exp(-(((xi - xp)/sigma)^2d0 + ((zi-zp)/sigma)^2d0))*2d0*(zi - zp)/sigma
      endforeach
    endforeach
    
    Um   = Um    + Um00*ampl*hz
    pres = pres  + gamm*p00*ampl*hz
    vx   = vx    + cs00*ampl*hz_vx
    vz   = vz    + cs00*ampl*hz_vz
    
    ; this calculates the new values of the 'densities' for momentum and energy:
    Upx  = Um*vx
    Upz  = Um*vz
    Ue   = pres/(gamm-1d0) + Um*(vx^2d0+vz^2d0)/2d0

  endif
end

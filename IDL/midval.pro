; ----------------------------------------------------------------------
; ROUTINE MIDVAL
;
;    PURPOSE: Calculate the averages of four neighbours grid points
;
;    INPUT ARGUMENTS: gg :matrix with arbitrary dimensions
;
;    OUTPUT: mid : a matrix with the mean of the 4 neiborg points
;
; ----------------------------------------------------------------------
function midval,gg
  return, (gg[0:-2,0:-2] + gg[1:-1,0:-2] + gg[0:-2,1:-1] + gg[1:-1,1:-1])/4d0
end

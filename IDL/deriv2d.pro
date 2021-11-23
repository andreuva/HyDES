; ----------------------------------------------------------------------
; ROUTINE DERIV2D
;
;    PURPOSE: Calculate the X and Z derivates for a matrix
;
;    INPUT ARGUMENTS: gg :matrix with arbitrary dimensions
;
;    OUTPUT: derx,derz : the derivates in either x or z directions 
;
; ----------------------------------------------------------------------
pro deriv2d,gg ,derx = derx, derz = derz

  ;call the commom blocks
  @grid_c
  
  ;compute the x or z derivates in each case or raise a error if the function is not properly call.
  if keyword_set(derx) then begin
    derx = ( (gg[1:-1,0:-2] - gg[0:-2,0:-2])/dx + (gg[1:-1,1:-1] - gg[0:-2,1:-1])/dx )/2d0
  endif
  
  if keyword_set(derz) then begin
    derz = ( (gg[0:-2,1:-1] - gg[0:-2,0:-2])/dz + (gg[1:-1,1:-1] - gg[1:-1,0:-2])/dz )/2d0
  endif 

  if (not keyword_set(derx)) and (not keyword_set(derz)) then begin
    print, 'deriv2d error : The derivates must be in the x or z direction (derx=.. or derz=.. as arguments)'
  endif
  
end

; ----------------------------------------------------------------------
; ROUTINE boundc
;
;    PURPOSE:  Calculate the boundary conditions of a variable
;    INPUT ARGUMENTS:  
;         typeofbc  : type of boundary conditions
;         site      : direcction in wich the boundary conditions are going to be set
;         varr      : variable to calculate the boundary conditions (will be changed)
;
;    OUTPUT:  The output is preformed via input arguments changing varr variable
;
;    NOTE: 'left' and 'right' below refer to the respective limits of the
;                 domain.
;          'periodic' is the type of boundary condition.
;          In the first practical, 'periodic' is the only type used, but
;          later on you may need to add further types of boundary conditions,
;          like, zero derivative, or linear extrapolation, or zero value, etc
;          etc.
; ----------------------------------------------------------------------

PRO boundc,typeofbc,site,varr

; compute the periodec boundary conditions either left or right
if typeofbc eq 'periodic' then begin
  
  if site eq 'left'  then begin
    varr[0,0:-1] = varr[-2,0:-1]
  endif else if site eq 'right' then begin
    varr[-1,0:-1] =  varr[1,0:-1]
  endif else if site eq 'down' then begin
    varr[0:-1,0] =  varr[0:-1,-2]
  endif else if site eq 'up' then begin
    varr[0:-1,-1] =  varr[0:-1,1]
  endif else begin
    print, 'boundc_2d.pro: Periodic boundary conditions should be left, right, up or down'
  endelse
  
endif else if typeofbc eq '0deriv' then begin
  
  if site eq 'left'  then begin
    varr[0,0:-1] = varr[1,0:-1]
  endif else if site eq 'right' then begin
    varr[-1,0:-1] =  varr[-2,0:-1]
  endif else if site eq 'down' then begin
    varr[0:-1,0] =  varr[0:-1,1]
  endif else if site eq 'up' then begin
    varr[0:-1,-1] =  varr[0:-1,-2]
  endif else begin
    print, 'boundc_2d.pro: 0deriv boundary conditions should be left, right, up or down'
  endelse
  
endif else begin
  print, 'boundc_2d : Type of boundary conditions not properly set'
endelse

end


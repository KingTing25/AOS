PRO calc_binary_pos, P_, tt, x1_, x2_

;print,a,e,inc,Omega,omgf,T0,Mtot,mu1,mu2
  a = P_(0)
  e = P_(1)
  inc = P_(2)
  Omega = P_(3)
  omgf = P_(4)
  T0 = P_(5)
  M1 = P_(6)
  M2 = P_(7)
  Mtot = M1+M2
  mu1 = M1/Mtot
  mu2 = M2/Mtot
  G = 6.67d-8
  Torb = sqrt(4.*!PI^2./(G*Mtot)*a^3.)
  n = 2.*!PI/Torb
  MM = n*(tt-T0)
  MM_p = MM mod (2.*!PI)
  MM = (MM_p+(2.*!PI)) mod (2.*!PI)
  if (MM gt !PI) then begin
    sgnMM = -1.0
  endif else begin
    sgnMM = 1.0
  endelse
  EE = MM+sgnMM*0.85*e
  eerror=EE-MM
  while (eerror ne 0) do begin
      EE1 = EE -(EE-e*sin(EE)-MM)/(1.-e*cos(EE))
      eerror=EE1-EE
      print, (EE1-EE)/EE
      EE = EE1
  endwhile
  rr = a*(1.-e*cos(EE))
  ff = 2.*atan(sqrt((1.+e)/(1.-e))*tan(EE/2.))

  
  if (finite(total(EE)) eq 0) then begin
      print, 'ERROR in calculating E'
      stop
  endif
  XX = rr*(cos(Omega)*cos(omgf+ff)-sin(Omega)*sin(omgf+ff)*cos(inc))
  YY = rr*(sin(Omega)*cos(omgf+ff)+cos(Omega)*sin(omgf+ff)*cos(inc))
  ZZ = rr*sin(omgf+ff)*sin(inc)
  
  x1_=[XX,YY,ZZ]*mu2
  x2_=[XX,YY,ZZ]*mu1*(-1.)
  return

END

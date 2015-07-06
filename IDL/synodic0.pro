PRO synodic0, Nt
  n=0
  tt=dblarr(Nt+1)
  Energy=dblarr(Nt+1)
  Errorx=dblarr(Nt+1)
  Errorv=dblarr(Nt+1)
  Cj=dblarr(Nt+1)
  M1=2.d33
  M2=2.d33
  G=6.67d-8
  a=5*1.5d13
  Torb=(4.*!PI^2./(G*(M1+M2)))^0.5 * a^(3.0/2.)
  r=a
  w=2*!PI/Torb
  e=0.
  i=0.
  Om=0.
  omgf=0.
  T0=0.
  P_=dblarr(8)
  P_(0)=a
  P_(1)=e
  P_(2)=i
  P_(3)=Om
  P_(4)=omgf
  P_(5)=T0
  P_(6)=M1
  P_(7)=M2
  Mtot = M1+M2
  mu1 = M1/Mtot
  mu2 = M2/Mtot
  r_scale = 0.1
  r2 = r_scale*a
  r2_ = dblarr(Nt+1)
  r2_(0) = r2
  calc_binary_pos, P_, tt(0), x1_, x2_
  v_scale = 0.9
  ;x_=[cos(60./!radeg),sin(!PI/3),0]*a
  x_=x2_ + r2*[-1,0,0]
  ;v_=1.01*[-sin(!PI/3),cos(!PI/3),0]*2*!PI*a/Torb
  ;x_=x2_+[-r2,0.,0.]
  v_=[0., v_scale*(-sqrt(G*M2/r2)) - sqrt(G*Mtot/a), 0.] ;prograde/retrograde
  ;v_=[0.,-sqrt(G*Mtot/(2.*a)), 0.] ;circumbinary
  
  ;L1
  ;x3_=[0d,0d,0d]
  ;L2
  ;x3_=[-x2temp-5.23805d13,0d,0d]
  ;L3
  ;x3_=[x1temp+5.23805d13,0d,0d]
  ;L4
  ;x_=[0d,(3^0.5)*a/2,0d]
  ;L5
  ;x3_=[0d,-(3^0.5)*a/2,0d]
  
  nn=2*!PI/Torb
  Norb=0.5
  dt=Norb*Torb/Nt
  r13=total(((x_-x1_)*(x_-x1_)))^0.5
  r23=total(((x_-x2_)*(x_-x2_)))^0.5
  ;L1
  ;v3_=[0d,0d,0d]
  ;v3_=v3_+0.5*v1_
  ;L2
  ;v3_=[0d,(G*(M1)/(2*r))^0.5*x3_(0)/x2temp,0d]
  ;L3
  ;v3_=[0d,(G*(M2)/(2*r))^0.5*x3_(0)/x1temp,0d]
  ;L4
  ;v_=[-(G*(M1)/(2*r))^0.5*3^0.5,0d,0d]
  ;L5
  ;v3_=[(G*(M1)/(2*r))^0.5*3^0.5,0d,0d]
  
  tt(0)=0
  Energy(0) = 0.5*total(abs(v_*v_))-G*M1/total(((x_-x1_)*(x_-x1_)))^0.5-G*M2/total(((x_-x2_)*(x_-x2_)))^0.5
  a2=0.2;this is not a comment
  a3=0.3
  a4=0.6
  a5=1.
  a6=0.875
  b21=0.2
  b31=0.075
  b32=0.225
  b41=0.3
  b42=-0.9
  b43=1.2
  b51=-11./54.
  b52=2.5
  b53=-70./27.
  b54=35./27.
  b61=1631./55296.
  b62=175./512.
  b63=575./13824.
  b64=44275./110592.
  b65=253./4096.
  c1=37./378.
  c2=0.
  c3=250./621.
  c4=125./594.
  c5=0.
  c6=512./1771.
  cs1=2825./27648.
  cs2=0.
  cs3=18575./48384.
  cs4=13525./55296.
  cs5=277./14336.
  cs6=0.25
  openw,1,'output.dat'
  
  eps=1d-6 ;error tolerance
  s=0.9
  fail=0
  win=0
  small=0
  big=0
  it=0.
  while (it lt Nt) do begin
    ;while ((total(x1_) gt -3.7218d13) and (it lt Nt)) do begin
    ;am=[[cos(n/360),-(sin(n/360)),0d],[sin(n/360),cos(n/360),0d],[0d,0d,1d]]
    ;bm=[[sin(n/360),cos(n/360),0d],[-(cos(n/360)),sin(n/360),0d],[0d,0d,1d]]
    h=dt
    calcfs, x_,v_,P_,tt(it),w,x_p,v_p
    xk1_n = dt*x_p
    vk1_n = dt*v_p
    
    calcfs, x_ + b21*xk1_n, v_ + b21*vk1_n,P_,tt(it)+a2*dt,w,x_p,v_p
    xk2_n = dt*x_p
    vk2_n = dt*v_p
    
    calcfs, x_ + b31*xk1_n + b32*xk2_n, v_ + b31*vk1_n + b32*vk2_n,P_,tt(it)+a3*dt,w,x_p,v_p
    xk3_n = dt*x_p
    vk3_n = dt*v_p
    
    calcfs, x_ + b41*xk1_n + b42*xk2_n + b43*xk3_n, v_ + b41*vk1_n + b42*vk2_n + b43*vk3_n,P_,tt(it)+a4*dt,w,x_p,v_p
    xk4_n = dt*x_p
    vk4_n = dt*v_p
    
    calcfs, x_ + b51*xk1_n + b52*xk2_n + b53*xk3_n + b54*xk4_n, v_ + b51*vk1_n + b52*vk2_n + b53*vk3_n + b54*vk4_n,P_,tt(it)+a5*dt,w,x_p,v_p
    xk5_n = dt*x_p
    vk5_n = dt*v_p
    
    calcfs, x_ + b61*xk1_n + b62*xk2_n + b63*xk3_n + b64*xk4_n + b65*xk5_n, v_ + b61*vk1_n + b62*vk2_n + b63*vk3_n + b64*vk4_n + b65*vk5_n,P_,tt(it)+a6*dt,w,x_p,v_p
    xk6_n = dt*x_p
    vk6_n = dt*v_p
    
    ;RK5
    x_n = c1*xk1_n + c2*xk2_n + c3*xk3_n + c4*xk4_n + c5*xk5_n + c6*xk6_n + x_
    v_n = c1*vk1_n + c2*vk2_n + c3*vk3_n + c4*vk4_n + c5*vk5_n + c6*vk6_n + v_
    
    ;RK4, embedded
    xs_n = cs1*xk1_n + cs2*xk2_n + cs3*xk3_n + cs4*xk4_n + cs5*xk5_n + cs6*xk6_n + x_
    vs_n = cs1*vk1_n + cs2*vk2_n + cs3*vk3_n + cs4*vk4_n + cs5*vk5_n + cs6*vk6_n + v_
    
    ;scales
    
    ;truncation error
    Errorx(it) = total(((x_n - xs_n)^2))^0.5
    Errorv(it) = total(((v_n - vs_n)^2))^0.5
    ;adaptive stepsize control
    errmax=0.0
    if errmax lt abs(Errorx(it)/sqrt(total(x_n^2))) then errmax=abs(Errorx(it)/sqrt(total(x_n^2)))
    if errmax lt abs(Errorv(it)/sqrt(total(v_n^2))) then errmax=abs(Errorv(it)/sqrt(total(v_n^2)))
    ;errmax=errmax/eps
    if errmax gt eps then begin ;failed step
      ;if h/10 gt s*h*errmax^(-0.25) then h = h*0.1 else
      h=s*h*(eps/errmax)^(0.25)
      if h eq 0. then begin
        print,'Stepsize of 0.'
        stop
      endif
      dt=h
      fail=fail+1
    endif else begin ;successful step
      win=win+1
      ;if errmax gt 1.89d-4 then begin ;next step will be small
      ;if errmax gt 1.89d then begin ;next step will be small
      h=s*h*(eps/errmax)^(0.2)
      ;small=small+1
      ;endif else begin ;next step will be big
      ;h=5.*h
      ;big=big+1
      ;endelse
      dt=h
      x_ = x_n
      v_ = v_n
      tt(it+1)=tt(it)+dt
      it=it+1
      printf,1,tt(it),x1_(0),x1_(1),x1_(2),x2_(0),x2_(1),x2_(2),x_(0),x_(1),x_(2)
      ;n=n+1
      calc_binary_pos, P_,tt(0),x1_,x2_
      r2=total(((x_-x2_)*(x_-x2_)))^0.5
      r1=total(((x_-x1_)*(x_-x1_)))^0.5
      r2_(it) = r2
      ;x_=matrix_multiply(am,x_)
      ;v_=matrix_multiply(am,v_)-n*matrix_multiply(bm,x_)
      ;Energy(it) = 0.5*total(v_*v_)-G*M1/total(((x_-x1_)*(x_-x1_)))^0.5-G*M2/total(((x_-x2_)*(x_-x2_)))^0.5
      Cj(it) = 2*(G*M1/r1 + G*M2/r2) + 2*nn*(x_(0)*v_(1)-x_(1)*v_(0))-v_(0)^2-v_(1)^2-v_(2)^2
    endelse
    ;print,errmax
    w=w+1
  endwhile
  close,1
  print,Cj(5),Cj(it)
  print,Cj(it)-Cj(5),(Cj(it)-Cj(5))/Cj(5)
  !P.charsize = 1.5
  window,1
  plot,tt/Torb,(Cj-mean(Cj))/mean(Cj),xrange=[tt(5),tt(it)]/Torb
  window,2
  plot,tt/Torb,deriv(tt),xrange=[0.,tt(it)/Torb]
  window,3
  plot,tt/Torb,r2_,xrange=[0,tt(it)/Torb]
  STOP
END

PRO calcfs, x_,v_,P_,tt,w,x_p,v_p
  G=6.67d-8
  calc_binary_pos,P_,tt(0),x1_,x2_
  Mtot=G*(P_(6)+P_(7))
  mu1 = P_(6)/Mtot
  mu2 = P_(7)/Mtot
  x_p=v_
  r23=total(((x_-x2_)*(x_-x2_)))^0.5
  r13=total(((x_-x1_)*(x_-x1_)))^0.5
  v_p=[2*w*v_(1) + w*w*x_(0) - (mu1*(x_(0)+mu2)/r13^3 + mu2*(x_(0)-mu1)/r23^3),-2*w*v_(0)+w*w*x_(1) - (mu1/r13^3+mu2/r23^3)*x_(1),-(mu1/r13^3+mu2/r23^3)*x_(2)] 
END
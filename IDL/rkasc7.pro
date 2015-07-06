PRO rkasc7
Nt=1000
tt=dblarr(Nt+1)
Energy=dblarr(Nt+1)
Errorx=dblarr(Nt+1)
Errorv=dblarr(Nt+1)
M1=2d33
M2=2d33
M3=2d30
G=6.67d-8
a=5*1.5d13
r=a
x2temp = M1*a/(M1+M2)
x1temp = a - x2temp

x1_=[x1temp,0d,0d]
x2_=[-x2temp,0d,0d]
x3_=x2_*0.9
;L1
;x3_=[0d,0d,0d]
;L2
;x3_=[-x2temp-5.23805d13,0d,0d]
;L3
;x3_=[x1temp+5.23805d13,0d,0d]
;L4
;x3_=[0d,(3^0.5)*a/2,0d]
;L5
;x3_=[0d,-(3^0.5)*a/2,0d]

Torb=(4.*!PI^2./(G*(M1+M2)))^0.5 * a^(3.0/2.)
Norb=0.5
dt=Norb*Torb/Nt
r13=total(((x3_-x1_)*(x3_-x1_)))^0.5
r23=total(((x3_-x2_)*(x3_-x2_)))^0.5

v1_=[0d,(G*(M2)/(2*r))^0.5,0d]
v2_=[0d,-(G*(M1)/(2*r))^0.5,0d]
;v3_=[0d,(G*M1/(-x2_(0)/10))^0.5,0d]
;L1
v3_=[0d,0d,0d]
;v3_=v3_+0.5*v1_
;L2
;v3_=[0d,(G*(M1)/(2*r))^0.5*x3_(0)/x2temp,0d]
;L3
;v3_=[0d,(G*(M2)/(2*r))^0.5*x3_(0)/x1temp,0d]
;L4
;v3_=[-(G*(M1)/(2*r))^0.5*3^0.5,0d,0d]
;L5
;v3_=[(G*(M1)/(2*r))^0.5*3^0.5,0d,0d]

tt(0)=0
Energy(0) = 0.5*M1*total(v1_*v1_)+0.5*M2*total(v2_*v2_)+0.5*M3*total(v3_*v3_)-G*M1*M2/total(((x2_-x1_)*(x2_-x1_)))^0.5-G*M1*M3/total(((x3_-x1_)*(x3_-x1_)))^0.5-G*M2*M3/total(((x3_-x2_)*(x3_-x2_)))^0.5
;a2=0.2
;a3=0.3
;a4=0.6
;a5=1.
;a6=0.875
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

eps=1d-10
s=0.9
fail=0
win=0
small=0
big=0
it=0
while (it lt Nt) do begin
;while ((total(x1_) gt -3.7218d13) and (it lt Nt)) do begin
  it=it+1
  h=dt
  calcf, x1_,v1_,x2_,v2_,x3_,v3_,x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
  x1k1_n = dt*x1_p
  x2k1_n = dt*x2_p
  x3k1_n = dt*x3_p
  v1k1_n = dt*v1_p
  v2k1_n = dt*v2_p
  v3k1_n = dt*v3_p
  
  calcf, x1_ + b21*x1k1_n, v1_ + b21*v1k1_n, x2_ + b21*x2k1_n, v2_ + b21*v2k1_n, x3_ + b21*x3k1_n, v3_ + b21*v3k1_n, x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
  x1k2_n = dt*x1_p
  x2k2_n = dt*x2_p
  x3k2_n = dt*x3_p
  v1k2_n = dt*v1_p
  v2k2_n = dt*v2_p
  v3k2_n = dt*v3_p
  
  calcf, x1_ + b31*x1k1_n + b32*x1k2_n, v1_+ b31*v1k1_n + b32*v1k2_n, x2_ + b31*x2k1_n + b32*x2k2_n, v2_ + b31*v2k1_n + b32*v2k2_n, x3_ + b31*x3k1_n + b32*x3k2_n, v3_ + b31*v3k1_n + b32*v3k2_n, x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
  x1k3_n = dt*x1_p
  x2k3_n = dt*x2_p
  x3k3_n = dt*x3_p
  v1k3_n = dt*v1_p
  v2k3_n = dt*v2_p
  v3k3_n = dt*v3_p
  
  calcf, x1_ + b41*x1k1_n + b42*x1k2_n + b43*x1k3_n, v1_+ b41*v1k1_n + b42*v1k2_n + b43*v1k3_n, x2_ + b41*x2k1_n + b42*x2k2_n + b43*x2k3_n, v2_ + b41*v2k1_n + b42*v2k2_n + b43*v2k3_n, x3_ + b41*x3k1_n + b42*x3k2_n + b43*x3k3_n, v3_ + b41*v3k1_n + b42*v3k2_n + b43*v3k3_n, x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
  x1k4_n = dt*x1_p
  x2k4_n = dt*x2_p
  x3k4_n = dt*x3_p
  v1k4_n = dt*v1_p
  v2k4_n = dt*v2_p
  v3k4_n = dt*v3_p
  
  calcf, x1_ + b51*x1k1_n + b52*x1k2_n + b53*x1k3_n + b54*x1k4_n, v1_+ b51*v1k1_n + b52*v1k2_n + b53*v1k3_n + b54*v1k4_n, x2_ + b51*x2k1_n + b52*x2k2_n + b53*x2k3_n + b54*x2k4_n, v2_ + b51*v2k1_n + b52*v2k2_n + b53*v2k3_n + b54*v2k4_n, x3_ + b51*x3k1_n + b52*x3k2_n + b53*x3k3_n + b54*x3k4_n, v3_ + b51*v3k1_n + b52*v3k2_n + b53*v3k3_n + b54*v3k4_n, x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
  x1k5_n = dt*x1_p
  x2k5_n = dt*x2_p
  x3k5_n = dt*x3_p
  v1k5_n = dt*v1_p
  v2k5_n = dt*v2_p
  v3k5_n = dt*v3_p
  
  calcf, x1_ + b61*x1k1_n + b62*x1k2_n + b63*x1k3_n + b64*x1k4_n + b65*x1k5_n, v1_+ b61*v1k1_n + b62*v1k2_n + b63*v1k3_n + b64*v1k4_n + b65*v1k5_n, x2_ + b61*x2k1_n + b62*x2k2_n + b63*x2k3_n + b64*x2k4_n + b65*x2k5_n, v2_ + b61*v2k1_n + b62*v2k2_n + b63*v2k3_n + b64*v2k4_n + b65*v2k5_n, x3_ + b61*x3k1_n + b62*x3k2_n + b63*x3k3_n + b64*x3k4_n + b65*x3k5_n, v3_ + b61*v3k1_n + b62*v3k2_n + b63*v3k3_n + b64*v3k4_n + b65*v3k5_n, x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
  x1k6_n = dt*x1_p
  x2k6_n = dt*x2_p
  x3k6_n = dt*x3_p
  v1k6_n = dt*v1_p
  v2k6_n = dt*v2_p
  v3k6_n = dt*v3_p
  
  ;RK5
  x1_n = c1*x1k1_n + c2*x1k2_n + c3*x1k3_n + c4*x1k4_n + c5*x1k5_n + c6*x1k6_n + x1_
  x2_n = c1*x2k1_n + c2*x2k2_n + c3*x2k3_n + c4*x2k4_n + c5*x2k5_n + c6*x2k6_n + x2_
  x3_n = c1*x3k1_n + c2*x3k2_n + c3*x3k3_n + c4*x3k4_n + c5*x3k5_n + c6*x3k6_n + x3_
  v1_n = c1*v1k1_n + c2*v1k2_n + c3*v1k3_n + c4*v1k4_n + c5*v1k5_n + c6*v1k6_n + v1_
  v2_n = c1*v2k1_n + c2*v2k2_n + c3*v2k3_n + c4*v2k4_n + c5*v2k5_n + c6*v2k6_n + v2_
  v3_n = c1*v3k1_n + c2*v3k2_n + c3*v3k3_n + c4*v3k4_n + c5*v3k5_n + c6*v3k6_n + v3_
  
  ;RK4, embedded
  x1s_n = cs1*x1k1_n + cs2*x1k2_n + cs3*x1k3_n + cs4*x1k4_n + cs5*x1k5_n + cs6*x1k6_n + x1_
  x2s_n = cs1*x2k1_n + cs2*x2k2_n + cs3*x2k3_n + cs4*x2k4_n + cs5*x2k5_n + cs6*x2k6_n + x2_
  x3s_n = cs1*x3k1_n + cs2*x3k2_n + cs3*x3k3_n + cs4*x3k4_n + cs5*x3k5_n + cs6*x3k6_n + x3_
  v1s_n = cs1*v1k1_n + cs2*v1k2_n + cs3*v1k3_n + cs4*v1k4_n + cs5*v1k5_n + cs6*v1k6_n + v1_
  v2s_n = cs1*v2k1_n + cs2*v2k2_n + cs3*v2k3_n + cs4*v2k4_n + cs5*v2k5_n + cs6*v2k6_n + v2_
  v3s_n = cs1*v3k1_n + cs2*v3k2_n + cs3*v3k3_n + cs4*v3k4_n + cs5*v3k5_n + cs6*v3k6_n + v3_
  
  ;scales
  
  ;truncation error
  Errorx(it) = total(((x1_n - x1s_n)^2+(x2_n - x2s_n)^2+(x3_n - x3s_n)^2))^0.5
  Errorv(it) = total(((v1_n - v1s_n)^2+(v2_n - v2s_n)^2+(v3_n - v3s_n)^2))^0.5
  ;adaptive stepsize control
  errmax=0.0
  errmax=max([abs(Errorx(it)/sqrt(total(x3_n^2))), abs(Errorv(it)/sqrt(total(v3_n^2)))])
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
    it=it-1
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
    x1_ = x1_n
    x2_ = x2_n
    x3_ = x3_n
    v1_ = v1_n
    v2_ = v2_n
    v3_ = v3_n
    printf,1,tt(it),x1_(0),x1_(1),x1_(2),x2_(0),x2_(1),x2_(2),x3_(0),x3_(1),x3_(2)
    tt(it)=tt(it-1)+dt
    Energy(it) = 0.5*M1*total(v1_*v1_)+0.5*M2*total(v2_*v2_)+0.5*M3*total(v3_*v3_)-G*M1*M2/total(((x2_-x1_)*(x2_-x1_)))^0.5-G*M1*M3/total(((x3_-x1_)*(x3_-x1_)))^0.5-G*M2*M3/total(((x3_-x2_)*(x3_-x2_)))^0.5
  endelse
  ;print,errmax
endwhile
close,1
print,Energy(0),Energy(it)
print,Energy(it)-Energy(0),(Energy(it)-Energy(0))/Energy(0)
window,1
plot,deriv(Energy),xrange=[0,it]
window,2
plot,deriv(tt),xrange=[0,it]
window,0
STOP 
END

PRO calcf, x1_k,v1_k,x2_k,v2_k,x3_k,v3_k,x1_p,v1_p,x2_p,v2_p,x3_p,v3_p,M1,M2,M3
G=6.67d-8
x1_p=v1_k
x2_p=v2_k
x3_p=v3_k
r12=total(((x2_k-x1_k)*(x2_k-x1_k)))^0.5
r23=total(((x3_k-x2_k)*(x3_k-x2_k)))^0.5
r13=total(((x3_k-x1_k)*(x3_k-x1_k)))^0.5
v1_p=G*M2/r12^3*(x2_k-x1_k)
v2_p=G*M1/r12^3*(x1_k-x2_k)
v3_p=G*M1/r13^3*(x1_k-x3_k) + G*M2/r23^3*(x2_k-x3_k)
END
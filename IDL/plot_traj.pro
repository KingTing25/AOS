PRO plot_traj,Nt,Torb,s
  window,0
  rdata = fltarr(10,Nt)
  tt = fltarr(Nt)
  x1 = fltarr(Nt)
  y1 = fltarr(Nt)
  z1 = fltarr(Nt)
  x2 = fltarr(Nt)
  y2 = fltarr(Nt)
  z2 = fltarr(Nt)
  x3 = fltarr(Nt)
  y3 = fltarr(Nt)
  z3 = fltarr(Nt)
  openr,1,'output.dat'
  readf,1,rdata
  close,1
  tt(*) = rdata(0,*)
  x1(*) = rdata(1,*)
  y1(*) = rdata(2,*)
  z1(*) = rdata(3,*)
  x2(*) = rdata(4,*)
  y2(*) = rdata(5,*)
  z2(*) = rdata(6,*)
  x3(*) = rdata(7,*)
  y3(*) = rdata(8,*)
  z3(*) = rdata(9,*)
  rr1 = sqrt(x1^2+y1^2+z1^2)
  rr2 = sqrt(x2^2+y2^2+z2^2)
  vx = deriv(tt,x3)
  vy = deriv(tt,y3)
  vz = deriv(tt,z3)
  rr = sqrt(x3^2+y3^2+z3^2)
  r1 = sqrt((x1-x3)^2+(y1-y3)^2+(z1-z3)^2)
  r2 = sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)
  rmax = max(rr)
  G = 6.67d-8
  M1 = 2d33
  M2 = 2d30
  nn = 2.*!PI/Torb
  Cj = 2.*(G*M1/r1 + G*M2/r2) + 2.*nn*(x3*vy-y3*vx)-vx^2.-vy^2.-vz^2.
  
  x1c = cos(nn*tt)*x1 + sin(nn*tt)*y1
  y1c = -sin(nn*tt)*x1 + cos(nn*tt)*y1
  z1c = z1
  x2c = cos(nn*tt)*x2 + sin(nn*tt)*y2
  y2c = -sin(nn*tt)*x2 + cos(nn*tt)*y2
  z2c = z2
  xc = cos(nn*tt)*x3 + sin(nn*tt)*y3
  yc = -sin(nn*tt)*x3 + cos(nn*tt)*y3
  zc = z3

  ;plot,[0,0],[0,0],xrange=[min(xc),max(xc)],yrange=[min(yc),max(yc)],$
  plot,[0,0],[0,0],xrange=[-1.2,1.2]*rmax,yrange=[-1.2,1.2]*rmax,$
  ;plot,[0,0],[0,0],xrange=[-1.2,-0.8]*rmax,yrange=[-0.2,0.2]*rmax,$
  ;plot,[0,0],[0,0],xrange=[-1.05,-0.95]*rmax,yrange=[-0.05,0.05]*rmax,$
    xstyle=1,ystyle=1,/isotropic

if s eq 0 then begin
  plots,x1,y1;,color=120
  plots,x2,y2,linestyle=1;,color=120
  plots,x3,y3,linestyle=2
endif else begin
  plots,x1c,y1c;,thick=4,color=120
  plots,x2c,y2c;,linestyle=1,thick=4,color=120
  plots,xc,yc
endelse
;  plots,xx/(aa/aa(0))/aa(0),yy/(aa/aa(0))/aa(0)

stop
END
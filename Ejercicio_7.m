clear; 
close all; 
zoom=5;
dr=0.1;
dt=2*pi/40;
r=3:dr:6;
t=0:dt:2*pi+dt;
[rr,tt]=meshgrid(r,t);
xx=rr.*cos(tt);
yy=rr.*sin(tt);
gg=2.*(rr-4./rr).*sin(tt);
Ur=tt;
Ut=tt;
Ux=1.*cos(tt);
Uy=1.*sin(tt);

figure;
axis([-zoom,zoom,-zoom,zoom])
plot(3*cos(t),3*sin(t),'k','linewidth',1)
hold on
plot(6*cos(t),6*sin(t),'k','linewidth',1)
hold on
% contour(xx,yy,gg,50);
quiver(xx,yy,Ux,Uy)
hold off
title('Campo electrico');
xlabel('x'); 
ylabel('y');
% plot(1:10)
% str = {'a = 3','b = 6'};
% text(2,7,str)
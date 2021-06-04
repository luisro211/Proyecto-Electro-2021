clear; 
close all; 
%-------------------------------------------------------------------------%
%                             INITIALIZATION                                  
%-------------------------------------------------------------------------%
% Constant 1/(4*pi*epsilon_0) = 9*10^9
k = 9*10^9;
% eps_r = Relative permittivity
eps_r = 1;
charge_order = 10^-9; % milli, micro, nano etc..
const = k*charge_order/eps_r;

% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
Nx = 50; % For 1 meter
Ny = 50; % For 1 meter
Nz = 50;
% n = Number of charges
n = 4;
% Electric fields Initialization
% E = Total electric field
E_f = zeros(Nx,Ny);
V   = zeros(Nx,Ny);     %voltage
% Ex = X-Component of Electric-Field
% Ey = Y-Component of Electric-Field
Ex = E_f;
Ey = E_f;
% ex = unit vector for x-component electric field
% ey = unit vector for y-component electric field
%ex = E_f;
%ey = E_f;
% r = distance between a selected point and the location of charge
r = E_f;
%r_square = E_f;
% Array of charges
% Q = All the 'n' charges are stored here
Q = [1,-1,1,-1];
% Array of locations
Xq = [-2,-2,2,2];

Yq = [2,-2,-2,2];
x_range = (1:Nx)-25;        %x coordinates for calculations
y_range = (1:Ny)-25;        %y coordinates for calculations
z_range = (1:Nz)-25;
[xMesh,yMesh]=meshgrid(x_range,y_range); %make 2 arrays with size (Nx,Ny)

%-------------------------------------------------------------------------%
%                      COMPUTATION OF ELECTRIC FIELDS
%-------------------------------------------------------------------------%

%  Repeat for all the 'n' charges
for k = 1:n
    r_square=(xMesh-Xq(k)).^2+(yMesh-Yq(k)).^2;
    r=sqrt(r_square);
    %extemp=const*Q(k)*(xMesh-Xq(k))./r;
    %eytemp=const*Q(k)*(yMesh-Yq(k))./r;
    %figure;
    %quiver(x_range,y_range,extemp,eytemp);
    %axis([-15 15 -15 15]); axis square;
    %tstr=sprintf('k=%d',k); title(tstr);
    Ex = Ex + const*Q(k)*4*(xMesh-Xq(k))./r.^3;      %x-component of field
    Ey = Ey + const*Q(k)*4*(yMesh-Yq(k))./r.^3;      %y-component of field
    V = V + const*Q(k)./r.^2;                      %voltage
end
E_f=sqrt(Ex.^2+Ey.^2);      %electric field magnitude
%-------------------------------------------------------------------------%
%                           PLOT THE RESULTS
%-------------------------------------------------------------------------%
figure;
% subplot(1,3,1);
contour_range = 0:1:20;
contour(x_range,y_range,E_f,contour_range,'linewidth',0.5);
title('E-field strength contours');
axis([-10 10 -10 10]);
axis square;
xlabel('y');
ylabel('z');

figure;
% subplot(1,3,2);
quiver(x_range,y_range,Ex,Ey);
title('Electric Field');
axis([-10 10 -10 10]);
axis square;
xlabel('y');
ylabel('z');

figure;
% subplot(1,3,3);
contour_range = -10:0.1:10;
contour(x_range,y_range,V,contour_range);
title('Voltage Contours');
axis([-10 10 -10 10]);
axis square;
xlabel('y');
ylabel('z');

figure;
surf(xMesh,yMesh,V);
title('Voltaje');
zlabel('Magnitud');
xlabel('Y'); 
ylabel('Z');

figure;
surf(x_range,y_range,E_f);
% contour(x_range,y_range,E_f,contour_range);
title('Campo electrico');
zlabel('Magnitud');
xlabel('Y'); 
ylabel('Z');
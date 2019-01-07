%%Parachute section model

clear
close all

%Radius at base
a = 6;
%Height including capped portion
b = 4;
%Hole or cap size
rmin = 1;
%Points plotted
numpoints = 30;
%Num of sections
numS = 16;

angmax = acos(rmin/a);
angle = linspace(0,angmax,numpoints);
cos(angle(1))
x = (2/numS)*pi*a*cos(angle);

fun = @(t) sqrt((a*sin(t)).^2 + (b*cos(t)).^2);

for i = 1:numpoints
    ang = angle(i);
    y(i) = integral(fun,0,ang);
end

dx = gradient(x, angle);
dy = gradient(y, angle);

magnitude = sqrt(dx.^2 + dy.^2);

xnew = x + .5*(dy./magnitude);
ynew = y - .5*(dx./magnitude);

xz = 0.*x;

plot(x,y,xnew,ynew,xz,y)
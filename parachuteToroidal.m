clear
close all

%% Parachute constants

Cd = 2.2;  %From Fruity Chute specs for comparable parachute and other research

%p = 1.225;  %Sea level
p = 1.00;   %For altitude of site ~4600 ft
p = .95;    %For altitude at main deployment ~6000 ft -> speed after main deployment

m = 40; % %lb
g = 9.8;    %m/s^2

%% Parachute size

v = 20;                 %Target velocity ft/s

v = v*12*2.54/100;      %m/s
m = m/2.2;              %kg
S = 2*m*g/(p*(v^2)*Cd); %Area of crossection at opening
R = sqrt(S/pi);         %Radius (m)
R = R*100/2.54          %inches

%% Size to velocity

R = 45;                     %Parachute radius in

R = R*2.54/100;             %m
S = pi*(R^2);               %Area of opening
v = sqrt(2*m*g/(p*S*Cd));   %Velocity m/s
v = v*100/(12*2.54)         %ft/s


%% Pattern inputs

%Total radius at base
A = 45; %18;
%Max height of parachute
b = 22.5; %9;
%Gap in center
C = 6.25; %2.5;
%Height diff between radii
h = sqrt(22.5^2 - 16.25^2);
%Points plotted
numpoints = 60;
%Num of sections
numS = 16;

%% Shroud lines

Ls = numS*(1.1*2*A/36 + .25) + numS*((sqrt(h^2 + C^2)/36 + .25)) + (numS/4)*(sqrt((1.1*2*A)^2 - A^2)/36  + .05) %yds total needed


%% Pattern
angmax = pi - asin(h/b)
a = (A - C)/(1 + abs(cos(angmax)));
c = A - 2*a;

angle = linspace(0,angmax,numpoints);
x = (1/numS)*pi*(a*(cos(angle) + 1) + c);

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

plot(x,y,xnew,ynew,zeros(size(y)),y,-x,y,-xnew,ynew)
grid on

angle(1)
x(1)
y(1)
%% Parachute image

for r = C:.2:A
    x = (-r):1:r;
    y = sqrt(r^2 - x.^2);
    X = [x,-x,x(1)];
    Y = [y,-y,y(1)];
    z = b*sqrt(1 - ((r - a - c)^2)/(a^2))*ones(size(X));
    
    figure(2)
    plot3(X,Y,z)
    hold on
end


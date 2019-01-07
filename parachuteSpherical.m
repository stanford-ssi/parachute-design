%%Parachute section model

clear
close all

%% Parachute constants

%Cd = 1.75;  %From Fruity Chute specs for comparable parachute and other research
Cd = 1.5;   %Conservative Cd used in design

%p = 1.225;  %Sea level
p = 1.00;   %For altitude of site ~4600 ft
p = .95;    %For altitude at main deployment ~6000 ft -> speed after main deployment

m = 36; % %lb
g = 9.8;    %m/s^2

%% Parachute size

v = 80;                 %Target velocity ft/s

v = v*12*2.54/100;      %m/s
m = m/2.2;              %kg
S = 2*m*g/(p*(v^2)*Cd); %Area of crossection at opening
R = sqrt(S/pi);         %Radius (m)
R = R*100/2.54          %inches

%% Size to velocity

R = 14;                     %Parachute radius in

R = R*2.54/100;             %m
S = pi*(R^2);               %Area of opening
v = sqrt(2*m*g/(p*S*Cd))   %Velocity m/s
v = v*100/(12*2.54)         %ft/s

%% Pattern design inputs

%Parachute radius inches
R = 70;
%Percent of area for vent
Per = 2;                    %Can be anywhere from 1 to 10
%Seam allowance
add = .5;

%Points plotted
numpoints = 50;
%Num of sections
numS = 24;

%% Shroud lines

Ls = numS*(1.1*2*R/36 + .25) %yds total needed

%% Parachute design

%Calculate vent area
So = 2*pi*(R^2);        
Av = So*Per/100;
rmin = sqrt(Av/(2*pi))

%Map segment
angmax = acos(rmin/R);
angle = linspace(0,angmax,numpoints);
x = (1/numS)*pi*R*cos(angle);
y = angle*R;

%Seam allowance
dx = gradient(x, angle);
dy = gradient(y, angle);

magnitude = sqrt(dx.^2 + dy.^2);

xnew = x + add*(dy./magnitude);
ynew = y - add*(dx./magnitude);

plot(x,y,xnew,ynew,zeros(size(y)),y,-x,y,-xnew,ynew)
grid on
xlabel("in")
ylabel("in")

%% Parachute image

N = 35;                      %Number of contour lines
step = (R - rmin)/(N - 1);

for r = rmin:step:R
    x = (-r):1:r;
    y = sqrt(r^2 - x.^2);
    X = [x,-x,x(1)];
    Y = [y,-y,y(1)];
    z = sqrt(R^2 - r^2)*ones(size(X));
    
    figure(2)
    plot3(X,Y,z)
    hold on
end


%% Graph Velocity for varying CD and m

clear

m = 5:5:60;
Cd = 1:.25:2.25;
R = 16;     %in
g = 9.8;    %m/s
p = 1.225;

for i = 1:numel(Cd)
   for j = 1:numel(m)
       v(j) = finalParaVSphere(R,m(j),g,p,Cd(i));
   end
   figure(3)
   plot(m,v)
   hold on
end

maxVal = 30;
max = maxVal*ones(size(m));

figure(3)
plot(m,max)
title("Weight versus Terminal Velocity")
xlabel("Weight (lb)")
ylabel("Velocity (ft/s)")
legend("Cd = " + string(Cd), "Max velocity")


function f = finalParaVSphere(R,m,g,p,Cd)
    %Returns velocity in ft/s
    
    R = R*2.54/100;         %m
    S = pi*(R^2);           %Area of opening
    m = m/2.2;              %kg
    v = sqrt(2*m*g/(p*S*Cd)); %Velocity m/s
    f = v*100/(12*2.54);
end







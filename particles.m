% particles.m
% Simulation of particles in a box. This is a more generic version of
% algorithms written by Rodney Metoyer and Chris Yoder for MAE 513.

% Begin by clearing the Matlab workspace. Also close all figures and clear
% the command window.
clear all; close all; clc;

% Simulation definition
startTime = 0.0;    % Start time in seconds
endTime = 30.0;     % End time in seconds
timeStep = 0.1;     % Output time step if desired
stateFunc = @state;
doYouWantMovie = true;
movieFile = 'test10.avi';
frameRate = 20;
speedReduction = 1.0;

% The vertices of an n-dimensional cube defines the space
height = 2.0;
width = 2.0;
space.box = [0 width width  0;...
             0 0     height height];
     
% Define the particles
particle.number = 10;   % Number of particles - must be an integer
particle.number = int32(particle.number); % Let's not take any chances. Note that int32 rounds, does not truncate
particle.radius = NaN(1,particle.number);
particle.mass = NaN(1,particle.number);
particle.spring = NaN(1,particle.number);
radius = 0.1; % For a homogenous radius distribution
mass = 0.1;
spring = 100.1;
for i=1:1:particle.number
    particle.radius(i) = radius;
    particle.mass(i) = mass;
    particle.spring(i) = spring;
end
% Change each radius individually if you want
particle.mass(5) = 0.1;

% Other environmental conditions
space.gravity = 1.0;

% Particle initial conditions
% Change initial conditions for each particle individually if you like
% These are just spaced evenly with the same speed.
% WARNING - current method assumes that you did not saturate your space!!!
xx0 = linspace(1.1*radius, width-1.1*radius,particle.number);
xy0 = ones(1,particle.number)*(height-1.1*radius);
xxd0 = 0.1*ones(1,particle.number);
xyd0 = zeros(1,particle.number);

% Put all initial conditions into one vector
for i = 1:1:particle.number    
    x0(4*(i - 1) + 1) = xx0(i);             % insert x conditions
    x0(4*(i - 1) + 2) = xy0(i);             % insert y conditions
    x0(4*(i - 1) + 3) = xxd0(i);            % insert xd conditions
    x0(4*(i - 1) + 4) = xyd0(i);            % insert yd conditions
end

%times = startTime:timeStep:endTime;

%t = 0; x = x0;
%xdot = stateFunc(t,x,space,particle);
times = linspace(startTime,endTime,endTime*frameRate);
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[time, states] = ode45(@(t,x)stateFunc(t,x,space,particle),times,x0);

% Break out for plotting and movie
for i = 1:1:particle.number    
    x(:,i) = states(:,(4*(i - 1) + 1));
    y(:,i) = states(:,(4*(i - 1) + 2));
    xd(:,i) = states(:,(4*(i - 1) + 3));
    yd(:,i) = states(:,(4*(i - 1) + 4));
end

% Make a movie of the motion
if doYouWantMovie
    figure;
    movegui(gcf);
    clear Mov; % Just in case
    n = 0;
    for i=1:1:length(x)
        for j = 1:1:particle.number
            plot(x(i,j),y(i,j),'ob','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b');
            hold on
        end
        hold off
        axis([0,width,0,height]); grid off;
        Mov(i) = getframe(gcf);
    end
    writerObj = VideoWriter(movieFile);
    writerObj.FrameRate = frameRate/speedReduction; writerObj.Quality = 100; % optional
    open(writerObj); writeVideo(writerObj,Mov); close(writerObj);
end
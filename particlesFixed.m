% particles.m
% Simulation of particles in a box. This is a more generic version of
% algorithms written by Rodney Metoyer and Chris Yoder for MAE 513.

% Begin by clearing the Matlab workspace. Also close all figures and clear
% the command window.
clear all; close all; clc;      % clear the command line, close figures, clear variables

% Simulation definition
startTime = 0.0;        % Start time in seconds
endTime = 30.0;         % End time in seconds
timeStep = 0.1;         % Output time step if desired NOTE: THIS DOES NOT CHANGE THE SOLVER TIMESTEP!!!
stateFunc = @stateFixed;
doYouWantMovie = true;  % true = make a movie file
movieFile = 'test1.avi';% name of the movie file
frameRate = 20;         % frame rate of the movie file
speedReduction = 1.0;   % reduce the frame rate by a constant value

% The vertices of an n-dimensional cube defines the space
height = 10.0;           % meters
width = 10.0;            % meters
space.box = [0 width width  0;...   % dimensions of the box
             0 0     height height];
     
% Define the particles
particle.number = 10;   % Number of particles - must be an integer
particle.number = int32(particle.number); % Let's not take any chances. Note that int32 rounds, does not truncate
particle.radius = NaN(1,particle.number); % preallocate vector of particle radii
particle.mass = NaN(1,particle.number);   % preallocate vector of particle masses
particle.spring = NaN(1,particle.number); % preallocate vector of particle spring constants
particle.damper = NaN(1, particle.number);% preallocate vector of particle damper constants
radius = 0.2;  % m, For a homogenous radius distribution
mass = 0.1;     % kg, mass of particles
spring = 100;   % N/m, spring constant
damper = 0.0;     % kg/s, damper constant
for i=1:1:particle.number           % for all particles,
    particle.radius(i) = radius;    % radius
    particle.mass(i) = mass;        % mass
    particle.spring(i) = spring;    % spring constant
    particle.damper(i) = damper;    % damper constant
end
% Change particle properties individually if you want
% particle.mass(5) = 0.1;

% Other environmental conditions
space.gravity = 9.81;           % m/s2, gravity

% Particle initial conditions
% Change initial conditions for each particle individually if you like
% These are just spaced evenly with the same speed.
% WARNING - current method assumes that you did not saturate your space!!!
xx0 = linspace(2*radius, width-2*radius,particle.number);   % linearlly space particles in x
xy0 = ones(1,particle.number)*(height-1.1*radius);              % place all on the same height
xxd0 = 0.1*ones(1,particle.number);                             % x initial velocoity
xyd0 = zeros(1,particle.number);                                % y initial velocity

% Put all initial conditions into one vector
% vector = [x1, y1, xd1, yd1, x2, y2, .... , xn, yn, xdn, ydn]
for i = 1:1:particle.number    
    x0(4*(i - 1) + 1) = xx0(i);             % insert x conditions
    x0(4*(i - 1) + 2) = xy0(i);             % insert y conditions
    x0(4*(i - 1) + 3) = xxd0(i);            % insert xd conditions
    x0(4*(i - 1) + 4) = xyd0(i);            % insert yd conditions
end

% Change particle initial conditions individually if you like
% If you want a single odd particle make oddParticle particle number that
% you want to be odd, otherwise make it zero
oddParticle = 1;
%x0(oddParticle) = 1.5*x0(oddParticle);
%x0(oddParticle+1) = 1.5*x0(oddParticle+1);
x0(3) = x0(oddParticle+2)*5;
x0(4) = x0(oddParticle+3)*50;
particle.radius(oddParticle) = 0.2;
x0(1) = 0.5*x0(1);      % specify the x position of the first particle
x0(2) = 0.5*x0(2);      % specify the y position of the first particle

% You can use the timestep directly if you want. I like to calculate one
% based on the framerate that I want.
%times = startTime:timeStep:endTime;
times = linspace(startTime,endTime,endTime*frameRate);      % linear vector for time
options = odeset('RelTol',1e-9,'AbsTol',1e-9); % Solution times can go up pretty quickly if you turn the tolerance too low.
[time, states] = ode23(@(t,x)stateFunc(t,x,space,particle),times,x0);   % ode solver command

% Break out for plotting and movie
for i = 1:1:particle.number    
    x(:,i) = states(:,(4*(i - 1) + 1));     % creates vector of x positions
    y(:,i) = states(:,(4*(i - 1) + 2));     % creates vector of y positions
    xd(:,i) = states(:,(4*(i - 1) + 3));    % creates vector of x velocities
    yd(:,i) = states(:,(4*(i - 1) + 4));    % creates vector of y velocities
end

% Make a movie of the motion
if doYouWantMovie
    figure;
    movegui(gcf);
    clear Mov; % Just in case
    n = 0;
    size = 10;
    for i=1:1:length(x)
        for j = 1:1:particle.number
            if j == oddParticle && oddParticle ~= 0
                color = 'r';
                size = 20;
            else
                color = 'b';
                size = 10;
            end
            plot(x(i,j),y(i,j),'ob','MarkerSize',size,'MarkerEdgeColor',color,'MarkerFaceColor',color);
            hold on
        end
        hold off                        % turn off the plot
        axis([0,width,0,height]);       % set axis 
        grid on;                        % turn on the grid
        title(['Time = ', num2str(time(i)), ' seconds']); % put current time in the title
        Mov(i) = getframe(gcf);         % get the frame and compile it into the movie file
    end
    writerObj = VideoWriter(movieFile); % write the movie to a file
    writerObj.FrameRate = frameRate/speedReduction; writerObj.Quality = 100; % optional
    open(writerObj); writeVideo(writerObj,Mov); close(writerObj);
end

% Move the movie to the bin
% TODO(Rodney) gitignore the bin
if exist('bin','dir') == 0
    mkdir('bin');
end
movefile(movieFile,['bin/' movieFile]);
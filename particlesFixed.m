% particles.m
% Simulation of particles in a box. This is a more generic version of
% algorithms written by Rodney Metoyer and Chris Yoder for MAE 513.

% Begin by clearing the Matlab workspace. Also close all figures and clear
% the command window.
clear all; 
close all; 
clc;      % clear the command line, close figures, clear variables

% Simulation definition
startTime = 0.0;        % Start time in seconds
endTime = 90.0;         % End time in seconds
timeStep = 0.25;         % Output time step if desired NOTE: THIS DOES NOT CHANGE THE SOLVER TIMESTEP!!!
stateFunc = @state2FixedC;
doYouWantMovie = true;  % true = make a movie file
movieFile = 'capture2.avi';% name of the movie file
frameRate = 20;         % frame rate of the movie file
speedReduction = 1.0;   % reduce the frame rate by a constant value

% The vertices of an n-dimensional cube defines the space
height = 20;           % meters
width = 20;            % meters
space.box = [0 width width  0;...   % dimensions of the box
             0 0     height height];
     
% Define the particles
particle.number = 5;   % Number of particles - must be an integer
particle.number = int32(particle.number); % Let's not take any chances. Note that int32 rounds, does not truncate
particle.radius = NaN(1,particle.number); % preallocate vector of particle radii
particle.mass = NaN(1,particle.number);   % preallocate vector of particle masses
particle.spring = NaN(1,particle.number); % preallocate vector of particle spring constants
particle.damper = NaN(1, particle.number);% preallocate vector of particle damper constants
particle.charge = NaN(1, particle.number); % preallocate vector of particle charges
particle.ke = 8987551787.3681764;       % N m2 / C2, Coulombs constant
radius = 8.775e-16;     % m, radius of a proton 
mass = 1.6726219e-27;     % kg, mass of particles
spring = 0.0001;   % N/m, spring constant
damper = 0;
charge = 1.6022e-19;    % C, charge of an electron
for i=1:1:particle.number           % for all particles,
    particle.radius(i) = radius;    % radius
    particle.mass(i) = mass;        % mass
    particle.spring(i) = spring;    % spring constant
    particle.damper(i) = damper;    % damper constant
    particle.charge(i) = -charge;   % elementary charge 
end
particle.time = endTime;
particle.charge(1) = 5*charge;
particle.ke = 8987551787.3681764;

% Change particle properties individually if you want
space.gravity = 0;

% Particle initial conditions
% Change initial conditions for each particle individually if you like
% These are just spaced evenly with the same speed.
% WARNING - current method assumes that you did not saturate your space!!!

% Initial conditions for orbiting charges (n = 2)
xx0 = [0.5, 0.25, 0.75, 0.5, 0.5]*width;
xy0 = [0.5, 0.5, 0.5, 0.85, 0.15]*height;
xxd0 = [0, 0, 0, 0.3625, -0.3625];
xyd0 = [0, 0.3625, -0.3625, 0, 0];
% xyd0 = [0, 0.5, -0.5, 0, 0];

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

% You can use the timestep directly if you want. I like to calculate one
% based on the framerate that I want.
times = startTime:timeStep:endTime;
options = odeset('RelTol',1e-9,'AbsTol',1e-9); % Solution times can go up pretty quickly if you turn the tolerance too low.
tic
[time, states] = ode23(@(t,x)stateFunc(t,x,space,particle),times,x0);   % ode solver command
toc

% Break out for plotting and movie
for i = 1:1:particle.number    
    x(:,i) = states(:,(4*(i - 1) + 1));     % creates vector of x positions
    y(:,i) = states(:,(4*(i - 1) + 2));     % creates vector of y positions
    xd(:,i) = states(:,(4*(i - 1) + 3));    % creates vector of x velocities
    yd(:,i) = states(:,(4*(i - 1) + 4));    % creates vector of y velocities
end

% Make a movie of the motion
if doYouWantMovie
    f1 = figure;
    axis([0,width,0,height]);       % set axis 
    set(f1, 'Units', 'points');     % sets units to points
    h1 = get(f1, 'Position');       % gets position in units of points
    get(f1, 'Units');
    movegui(gcf);
    clear Mov; % Just in case
    n = 0;
    size = 10;
    for i=1:1:length(time)
        for j = 1:1:particle.number
            if j == oddParticle && oddParticle ~= 0
                color = 'r';
                %size = particle.radius(oddParticle)*h1(3)/width;
                size = 12;
            else
                color = 'b';
                %size = particle.radius(j)*h1(3)/width;
                size = 8;
            end
            plot(x(i,j),y(i,j),'ob','MarkerSize',size,'MarkerEdgeColor',color,'MarkerFaceColor',color);
            hold on
        end
        plot([space.box(1, :), space.box(1, 1)], [space.box(2, :),space.box(2, 1)], 'r') % plots the box
        axis([0,width,0,height]);       % set axis 
        axis equal
        grid on;                        % turn on the grid
        plot(states(1:i, 5), states(1:i, 6))
        plot(states(1:i, 9), states(1:i, 10))
        hold off
        %title(['Time = ', num2str(time(i)), ' seconds']); % put current time in the title
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
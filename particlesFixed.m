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
endTime = 20.0;         % End time in seconds
timeStep = 0.1;         % Output time step if desired NOTE: THIS DOES NOT CHANGE THE SOLVER TIMESTEP!!!
stateFunc = @state2FixedC;
%c1 = 0.3247; c2 = -3.1269; c3 = 10.372; % coefficients for stateFixed
c1 = 0.1999; c2 = -3.1515; c3 = 13.751; % coefficients for state2Fixed

doYouWantMovie = true;  % true = make a movie file
movieFile = 'capture1.avi';% name of the movie file
frameRate = 20;         % frame rate of the movie file
speedReduction = 1.0;   % reduce the frame rate by a constant value

% The vertices of an n-dimensional cube defines the space
height = 10;           % meters
width = 10;            % meters
space.box = [0 width width  0;...   % dimensions of the box
             0 0     height height];
     
% Define the particles
particle.number = 3;   % Number of particles - must be an integer
particle.number = int32(particle.number); % Let's not take any chances. Note that int32 rounds, does not truncate
particle.radius = NaN(1,particle.number); % preallocate vector of particle radii
particle.mass = NaN(1,particle.number);   % preallocate vector of particle masses
particle.spring = NaN(1,particle.number); % preallocate vector of particle spring constants
particle.damper = NaN(1, particle.number);% preallocate vector of particle damper constants
particle.charge = NaN(1, particle.number); % preallocate vector of particle charges
particle.ke = 8987551787.3681764;       % N m2 / C2, Coulombs constant
% radius = 0.05;  % m, For a homogenous radius distribution
radius = 8.775e-16;     % m, radius of a proton 
mass = 1.6726219e-27;     % kg, mass of particles
spring = 0.0001;   % N/m, spring constant
% damper = 2.5;     % kg/s, damper constant
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
particle.charge(1) = -2*charge;
particle.charge(2) = 8*charge;
particle.charge(3) = -2*charge;

% Change particle properties individually if you want
% particle.mass(5) = 0.1;
% Other environmental conditions
% space.gravity = 9.81;           % m/s2, gravity
space.gravity = 0;

% Particle initial conditions
% Change initial conditions for each particle individually if you like
% These are just spaced evenly with the same speed.
% WARNING - current method assumes that you did not saturate your space!!!

% xx0 = linspace(2*radius, width-2*radius,particle.number);   % linearlly space particles in x
% xy0 = ones(1,particle.number)*(height-1.1*radius);              % place all on the same height
xx0 = [0.25, 0.5, 0.75]*width;
xy0 = [0.5, 0.5, 0.5]*height;
% xx0(1:4) = [0.25, 0.75, 0.75, 0.25]*width;
% xy0(1:4) = [0.75, 0.75, 0.25, 0.25]*height;
% xx0(10:18) = 0.07*width*(1:9);
% xx0(19:27) = 0.07*width*(1:9);
% xx0(28:36) = 0.07*width*(1:9);
% xx0(37:45) = 0.07*width*(1:9);
% xx0(46:54) = 0.07*width*(1:9);
% xx0(55:63) = 0.07*width*(1:9);
% xx0(64:72) = 0.07*width*(1:9);
% xx0(73:81) = 0.07*width*(1:9);

% xy0(10:18) = 2*0.075*height;
% xy0(19:27) = 3*0.075*height;
% xy0(28:36) = 4*0.075*height;
% xy0(37:45) = 5*0.075*height;
% xy0(46:54) = 6*0.075*height;
% xy0(55:63) = 7*0.075*height;
% xy0(64:72) = 8*0.075*height;
% xy0(73:81) = 9*0.075*height;

xxd0 = zeros(1,particle.number);                             % x initial velocoity
xyd0 = zeros(1,particle.number);                                % y initial velocity

figure
plot(xx0, xy0, 'ob')
axis([0,width,0,height]);       % set axis
grid on

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
oddParticle = 0;
% particle.radius(oddParticle) = radius*3;
% particle.mass(oddParticle) = mass*3;
% x0(oddParticle) = 1.5*x0(oddParticle);
% x0(oddParticle+1) = 1.5*x0(oddParticle+1);
% x0(3) = x0(oddParticle+2)*5;
% x0(4) = x0(oddParticle+3)*50;
% particle.radius(oddParticle) = 0.2;
% x0(1) = 0.5*x0(1);      % specify the x position of the first particle
% x0(2) = 0.5*x0(2);      % specify the y position of the first particle

% You can use the timestep directly if you want. I like to calculate one
% based on the framerate that I want.
%times = startTime:timeStep:endTime;
tc = c1*(particle.number^2) + c2*particle.number + c3;
disp(['Estimated run time is ', num2str(tc), ' seconds']);
times = linspace(startTime,endTime,endTime*frameRate);      % linear vector for time
options = odeset('RelTol',1e-9,'AbsTol',1e-9); % Solution times can go up pretty quickly if you turn the tolerance too low.
tic
[time, states] = ode45(@(t,x)stateFunc(t,x,space,particle),times,x0);   % ode solver command
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
%     axis equal
    axis([0,width,0,height]);       % set axis 
    
    %h2 = get(f1, 'Position')       % gets position in units of points
    %disp('Pixels?');
    set(f1, 'Units', 'points');     % sets units to points
    h1 = get(f1, 'Position');       % gets position in units of points
    get(f1, 'Units');
    %rat = h1(4)/height;             % sets ratio of figure units to points
    %disp('Points?');
    %set(f1, 'Units', 'pixels');     % sets units back to pixels
    movegui(gcf);
    clear Mov; % Just in case
    n = 0;
    size = 10;
    for i=1:1:length(x)
        for j = 1:1:particle.number
            if j == oddParticle && oddParticle ~= 0
                color = 'r';
                %size = particle.radius(oddParticle)*h1(3)/width;
                size = 10;
            else
                color = 'b';
                %size = particle.radius(j)*h1(3)/width;
                size = 10;
            end
            plot(x(i,j),y(i,j),'ob','MarkerSize',size,'MarkerEdgeColor',color,'MarkerFaceColor',color);
            hold on
        end
        plot([space.box(1, :), space.box(1, 1)], [space.box(2, :), space.box(2, 1)], 'r')
        hold off                        % turn off the plot
%         axis equal
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
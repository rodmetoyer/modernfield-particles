function [xdot] = state(t,x,space,particle)
t                   % display the current time in the command line 
% init xdot
xdot = NaN(length(x),1);    % initialize the state derivative vector

% Particle independent properties
box = space.box;            % extract box size
g = space.gravity;          % extract box gravity


% Find the force on each particle
for i=1:1:particle.number   % for all the particles
    
    % Reset
    Fx = 0.0;               % set net force in x to zero
    Fy = 0.0;               % set net force in y to zero
    
    % Particle dependent properties
    radiusi = particle.radius(i);   % set the radius for the particle
    springi = particle.spring(i);   % set the spring for the particle
    massi = particle.mass(i);       % set the mass for the particle
    damperi = particle.damper(i);   % set the damper for the particle
    
    % Gravity
    Fy = -massi*g;                  % calculate gravity 
    
    % Walls
    if x(4*(i - 1) + 1) < radiusi + box(1,1)        % if x < radius + left wall, calculate force
        Fx = Fx + (radiusi - x(4*(i - 1) + 1))*springi + (-x(4*(i - 1) + 3))*damperi;
    elseif x(4*(i - 1) + 1) > box(1,2) - radiusi    % if x > right wall - radius, calculate force
        Fx = Fx - (radiusi - (box(1,2) - x(4*(i - 1) + 1)))*springi - (x(4*(i - 1) + 3))*damperi;
    end
    
    if x(4*(i - 1) + 2) < radiusi + box(2,1)        % if y < radius + bottom wall, calculate force
        Fy = Fy + (radiusi-x(4*(i - 1) + 2))*springi + (-x(4*(i - 1) + 4))*damperi;
    elseif x(4*(i - 1) + 2) > box(2,4) - radiusi    % if y > top wall - radius, calculate force
        Fy = Fy - (radiusi - (box(2,4) - x(4*(i - 1) + 2)))*springi - x(4*(i - 1) + 4)*damperi;
    end
    
    % Now do interparticle forces
    % Note that this is probably the worst way to do this, but it works
    for j=1:1:particle.number        
        % Don't let particles collide with themselves
        if i == j
            continue
        end
        
        % Get the radius of particle j
        radiusj = particle.radius(j);       % pull radius of the second particle
        
        % Get the distance between i and j
        distx = x(4*(j - 1) + 1) - x(4*(i - 1) + 1);
        disty = x(4*(j - 1) + 2) - x(4*(i - 1) + 2);
        distance = sqrt(distx^2 + disty^2);
        
        % Get the time rate of change of the distance between i and j
        dxd = x(4*(j - 1) + 3) - x(4*(i - 1) + 3);
        dyd = x(4*(j - 1) + 4) - x(4*(i - 1) + 4);
        distd = -(distx*dxd + disty*dyd)/distance; % derivative of R1 + R2 - distance
        
        % See if particle j is in contact with particle i. If it is,
        % calculate the force on i (force on j will be calculated when it
        % becomes i -> I told you this was a horrible way to do this).
        if distance < radiusi+radiusj
            % Force magnitude -> springs in series
            springj = particle.spring(j);
            damperj = particle.damper(j);
            Fmag = (radiusi+radiusj-distance)*(springi*springj)/(springi+springj);
            Fdmp = (-distd)*(damperi*damperj)/(damperi+damperj);
            Fx = Fx - Fmag*distx/distance + Fdmp*distx/distance; % Minus because it is force on i
            Fy = Fy - Fmag*disty/distance + Fdmp*distx/distance;
        end
    end
    
    % Put it all together
    xdot(4*(i - 1) + 1) = x(4*(i - 1) + 3);
    xdot(4*(i - 1) + 2) = x(4*(i - 1) + 4);
    xdot(4*(i - 1) + 3) = 1/massi*Fx;
    xdot(4*(i - 1) + 4) = 1/massi*Fy;   
    
    % debugging - to remove
    %disp(i);
    %xdot = x;
end
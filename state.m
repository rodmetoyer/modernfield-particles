function [xdot] = state(t,x,space,particle)

% init xdot
xdot = NaN(length(x),1);

% Particle independent properties
box = space.box;
g = space.gravity;

% Find the force on each particle
for i=1:1:particle.number
    
    % Reset
    Fx = 0.0;
    Fy = 0.0;
    
    % Particle dependent properties
    radiusi = particle.radius(i);
    springi = particle.spring(i);
    massi = particle.mass(i);
    
    % Gravity
    Fy = -massi*g;
    
    % Walls
    if x(4*(i - 1) + 1) < radiusi + box(1,1)
        Fx = Fx + (radiusi - x(4*(i - 1) + 1))*springi;
    elseif x(4*(i - 1) + 1) > box(1,2) - radiusi
        Fx = Fx - (radiusi - (box(1,2) - x(4*(i - 1) + 1)))*springi;
    end
    
    if x(4*(i - 1) + 2) < radiusi + box(2,1)
        Fy = Fy + (radiusi-x(4*(i - 1) + 2))*springi;
    elseif x(4*(i - 1) + 2) > box(2,4) - radiusi
        Fy = Fy - (radiusi - (box(2,4) - x(4*(i - 1) + 2)))*springi;
    end
    
    % Now do interparticle forces
    % Note that this is probably the worst way to do this, but it works
    for j=1:1:particle.number        
        % Don't let particles collide with themselves
        if i == j
            continue
        end
        
        % Get the radius of particle j
        radiusj = particle.radius(j);
        
        % Get the distance between i and j
        distx = x(4*(j - 1) + 1) - x(4*(i - 1) + 1);
        disty = x(4*(j - 1) + 2) - x(4*(i - 1) + 2);
        distance = sqrt(distx^2 + disty^2);
        
        % See if particle j is in contact with particle i. If it is,
        % calculate the force on i (force on j will be calculated when it
        % becomes i -> I told you this was a horrible way to do this).
        if distance < radiusi+radiusj
            % Force magnitude -> springs in series
            springj = particle.spring(j);
            Fmag = (radiusi+radiusj-distance)*(springi*springj)/(springi+springj);
            Fx = Fx - Fmag*distx/distance; % Minus because it is force on i
            Fy = Fy - Fmag*disty/distance;
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
function [xdot] = state(t,x,space,particle)

% init xdot
xdot = NaN(length(x),1);

% Find the force on each particle
for i=1:1:particle.number
    
    % Reset
    Fx = 0.0;
    Fy = 0.0;
    
    % Breakout to save typing
    radius = particle.radius(i);
    spring = particle.spring(i);
    mass = particle.mass(i);
    box = space.box;
    g = space.gravity;
    
    % Gravity
    Fy = -mass*g;
    
    % Walls
    if x(4*(i - 1) + 1) < radius + box(1,1)
        Fx = Fx + (radius - x(4*(i - 1) + 1))*spring;
    elseif x(4*(i - 1) + 1) > box(1,2) - radius
        Fx = Fx - (radius - (box(1,2) - x(4*(i - 1) + 1)))*spring;
    elseif x(4*(i - 1) + 2) < radius + box(2,1)
        Fy = Fy + (radius-x(4*(i - 1) + 2))*spring;
    elseif x(4*(i - 1) + 2) > box(2,4) - radius
        Fy = Fy - (radius - (box(2,4) - x(4*(i - 1) + 2)))*spring;
    end
    
    % Now do interparticle forces
    
    % Put it all together
    xdot(4*(i - 1) + 1) = x(4*(i - 1) + 3);
    xdot(4*(i - 1) + 2) = x(4*(i - 1) + 4);
    xdot(4*(i - 1) + 3) = 1/mass*Fx;
    xdot(4*(i - 1) + 4) = 1/mass*Fy;   
    
    % debugging - to remove
    %disp(i);
    %xdot = x;
end
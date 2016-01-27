function [xdot] = state2FixedC(t,x,space,particle)
% Fixed = fixed the issue with contact between particles  and the wall
% FixedC = adding Coulomb forces to the model (distance forces)

% init xdot
xdot = NaN(length(x),1);

% Particle independent properties
box = space.box;
g = space.gravity;

% Initalize matrices
FX = zeros(particle.number, particle.number);
FY = zeros(particle.number, particle.number);

% Find the force on each particle
for i=1:1:particle.number
    Fx = 0;
    Fy = 0;
    radiusi = particle.radius(i);
    springi = particle.spring(i);
    massi = particle.mass(i);
    damperi = particle.damper(i);
    chargei = particle.charge(i);
    
    % Walls
    if x(4*(i - 1) + 1) < radiusi + box(1,1) % if x < radius + left wall
        Fx = (box(1, 1) + radiusi - x(4*(i - 1) + 1))*springi + (-x(4*(i - 1) + 3))*damperi;
    elseif x(4*(i - 1) + 1) > box(1,2) - radiusi % if x > right wall - radius
        Fx = - (radiusi - (box(1,2) - x(4*(i - 1) + 1)))*springi - (x(4*(i - 1) + 3))*damperi;
    end
    
    if x(4*(i - 1) + 2) < radiusi + box(2,1) % if y < radius + bottom wall
        Fy = (box(2, 1) + radiusi - x(4*(i - 1) + 2))*springi + (-x(4*(i - 1) + 4))*damperi;
    elseif x(4*(i - 1) + 2) > box(2,4) - radiusi % if y > top wall - radius
        Fy = - (radiusi - (box(2,4) - x(4*(i - 1) + 2)))*springi - x(4*(i - 1) + 4)*damperi;
    end

% Gravity
Fy = -massi*g + Fy;

% wall forces
FX(i, i) = Fx;
FY(i, i) = Fy;
end

% MASS FORCES
for i=1:1:particle.number - 1   
    fx = 0;
    fy = 0;
    
    % Particle dependent properties
    radiusi = particle.radius(i);   % extract the radius of the particle
    springi = particle.spring(i);
    damperi = particle.damper(i);
    chargei = particle.charge(i);
    
    % Now do interparticle forces
    % Note that this is probably the second to worst way to do this, but it works
    for j=i+1:1:particle.number        
       
        % Get the radius of particle j
        radiusj = particle.radius(j);
        
        % Get the distance between i and j
        distx = x(4*(j - 1) + 1) - x(4*(i - 1) + 1);
        disty = x(4*(j - 1) + 2) - x(4*(i - 1) + 2);
        distance = sqrt(distx^2 + disty^2);
        
        % Get the time rate of change of the distance between i and j
        dxd = x(4*(j - 1) + 3) - x(4*(i - 1) + 3);
        dyd = x(4*(j - 1) + 4) - x(4*(i - 1) + 4);
        distd = -(distx*dxd + disty*dyd)/distance; % derivative of R1 + R2 - distance

        % Particle forces
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mechanical forces (contact for rigid body motion)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % See if particle j is in contact with particle i. If it is,
        % calculate the force on i (force on j will be calculated when it
        % becomes i -> I told you this was a horrible way to do this).
        if distance < radiusi+radiusj
            % Force magnitude -> springs in series
            springj = particle.spring(j);
            damperj = particle.damper(j);
            Fmag = (radiusi+radiusj-distance)*(springi*springj)/(springi+springj);
            if damperi + damperj < 1.0e-6
                Fdmp = 0.0;
            else
                Fdmp = (-distd)*(damperi*damperj)/(damperi+damperj);
            end
            fx = - Fmag*distx/distance + Fdmp*distx/distance; % Minus because it is force on i
%             FX = Fmag*distx/distance - Fdmp*distx/distance; % Minus because it is force on i
            fy = - Fmag*disty/distance + Fdmp*distx/distance;
%             FY = Fmag*disty/distance - Fdmp*distx/distance;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Electrical forces (Coulomb force)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get parameters
        ke = particle.ke;               % coulomb's constant
        chargej = particle.charge(j);   % charge on the particle
        
        % Calculate the force between particles
        Fe = ke*chargei*chargej/(distance^2);
        fx = fx - Fe*distx/distance;
        fy = fy - Fe*disty/distance;

        
        % Sum all forces on the particle
        FX(i, j) = -fx;      % force on i due to j
        FX(j, i) = fx;     % force on j due to i
        FY(i, j) = -fy;      % force on i due to j
        FY(j, i) = fy;     % force on j due to i
    end
end

for k = 1:1:particle.number   
    % Put it all together
    massi = particle.mass(k);
    xdot(4*(k - 1) + 1) = x(4*(k - 1) + 3);
    xdot(4*(k - 1) + 2) = x(4*(k - 1) + 4);
    xdot(4*(k - 1) + 3) = (1/massi)*sum(FX(k, :));
    xdot(4*(k - 1) + 4) = (1/massi)*sum(FY(k, :));   
end
end
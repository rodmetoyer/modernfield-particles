function [xdot] = stateCoulomb(t,x,space,particle)
% init xdot
xdot = NaN(length(x),1);

% Particle independent properties
box = space.box;
g = space.gravity;

% Initalize matrices
FX = zeros(particle.number, particle.number);
FY = zeros(particle.number, particle.number);

% Wall forces
for i=1:1:particle.number
    Fx = 0;
    Fy = 0;
    radiusi = particle.radius(i);
    springi = particle.spring(i);
    massi = particle.mass(i);
    damperi = particle.damper(i);
    
    % Walls
    if x(4*(i - 1) + 1) < radiusi + box(1,1) % if x < radius + left wall
        Fx = (radiusi - x(4*(i - 1) + 1))*springi + (-x(4*(i - 1) + 3))*damperi;
    elseif x(4*(i - 1) + 1) > box(1,2) - radiusi % if x > right wall - radius
        Fx = - (radiusi - (box(1,2) - x(4*(i - 1) + 1)))*springi - (x(4*(i - 1) + 3))*damperi;
    end
    
    if x(4*(i - 1) + 2) < radiusi + box(2,1) % if y < radius + bottom wall
        Fy = (radiusi-x(4*(i - 1) + 2))*springi + (-x(4*(i - 1) + 4))*damperi;
    elseif x(4*(i - 1) + 2) > box(2,4) - radiusi % if y > top wall - radius
        Fy = - (radiusi - (box(2,4) - x(4*(i - 1) + 2)))*springi - x(4*(i - 1) + 4)*damperi;
    end

    % Gravity
    Fy = -massi*g + Fy;
    
    % wall forces
    FX(i, i) = Fx;
    FY(i, i) = Fy;
end

% Interparticle FORCES
for i=1:1:particle.number - 1   
    
    ke = particle.ke;               % coulomb's constant
    
    % Particle dependent properties
    radiusi = particle.radius(i);
    springi = particle.spring(i);
    damperi = particle.damper(i);
    chargei = particle.charge(i);   
        
    for j=i+1:1:particle.number
        
        
        % Get the radius of particle j
        radiusj = particle.radius(j);

        % Get the distance between i and j
        distx = x(4*(j - 1) + 1) - x(4*(i - 1) + 1);
        disty = x(4*(j - 1) + 2) - x(4*(i - 1) + 2);
        distance = sqrt(distx^2 + disty^2);
        %TODO(Rodney) Need to build in some protection for coincidence when
        %when collisions are off.
        
        Fmag = 0.0;
        Fdmp = 0.0;
        
        if particle.collisions

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
                if damperi + damperj < 1.0e-6
                    Fdmp = 0.0;
                else
                    Fdmp = (-distd)*(damperi*damperj)/(damperi+damperj);
                end
            end
        end

        % Electrical forces (Coulomb force)        
        chargej = particle.charge(j);   % charge on the particle
        
        % Calculate the force between particles
        % fx is the Fx force in FX(i, j) which is the force on ith
        % particle with respect to the jth particle (ie the force
        % vector starts at j and ends at i)
        Fe = ke*chargei*chargej/(distance^2);
        fx = Fe*distx/distance;
        fy = Fe*disty/distance;
        
        FX(i, j) = - Fmag*distx/distance + Fdmp*distx/distance - fx; 
        FX(j, i) = Fmag*distx/distance - Fdmp*distx/distance + fx; 
        FY(i, j) = - Fmag*disty/distance + Fdmp*distx/distance - fy;
        FY(j, i) = Fmag*disty/distance - Fdmp*distx/distance + fy;        
    end
end

for k = 1:1:particle.number
    mass = particle.mass(k);
    % Put it all together
    xdot(4*(k - 1) + 1) = x(4*(k - 1) + 3);
    xdot(4*(k - 1) + 2) = x(4*(k - 1) + 4);
    xdot(4*(k - 1) + 3) = (1/mass)*sum(FX(k, :));
    xdot(4*(k - 1) + 4) = (1/mass)*sum(FY(k, :));   
end
end
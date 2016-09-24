clear

# initialize the code for a two body problem

load("finalData1.dat");

np = 10;  # number of objects
x = a(1:np,1);
y = a(1:np,2);
vx = a(1:np,3);
vy = a(1:np,4);
m = a(1:np,5);

# set the final times
tf = 1;
dt = 0.001;
nsteps = tf/dt;
ee = 0.05;  # softening length used to prevent close approach problems

# loop over the time range
for k = 1:nsteps;
  
  # zero the forces 
  fx = zeros(np,1);
  fy = zeros(np,1);
  pe = zeros(np,1);
  ke = zeros(np,1);
  
  # for each source particle
  for i = 1:np-1
    
    # loop over the other particles
    for j = i+1:np
      
      # calculate the displacements
      dx = x(i) - x(j);
      dy = y(i) - y(j);
      
      # calculate the distances
      dr2 = dx^2 + dy^2 + ee^2;
      dr = sqrt(dr2);
      
      # find the total force
      df = m(i) * m(j) / dr2;
      
      # resolve the force into components
      dfx = df * dx/dr;
      dfy = df * dy/dr;
      
      # update the target and source
      fx(i) = fx(i) - dfx;
      fy(i) = fy(i) - dfy; 
      fx(j) = fx(j) + dfx;
      fy(j) = fy(j) + dfy;
      
      # update potential energy
      pe(i) = pe(i) - m(i) *m(j) / dr /2;
      pe(j) = pe(j) - m(i) *m(j) / dr /2;
    end
  end
  
  # now move the particles
  for i = 1:np
    
    # calculate kinetic energy
    ke(i) = 0.5 * m(i) * (vx(i)^2 + vy(i)^2);
    
    # find the accelerations from the forces
    ax = fx(i) / m(i);
    ay = fy(i) / m(i);
    
    # find the change in velocity
    dvx = ax * dt;
    dvy = ay * dt;
    
    # update the velocities
    vx(i) = vx(i) + dvx;
    vy(i) = vy(i) + dvy;
    
    # find the change in position
    dx = vx(i) * dt ;
    dy = vy(i) * dt ;
    
    # update the position
    x(i) = x(i) + dx ;
    y(i) = y(i) + dy; 
  end
  
  # save the results
  pv1(k,1) = x(1);
  pv1(k,2) = y(1);
  pv1(k,3) = vx(1);
  pv1(k,4) = vy(1);
  pv1(k,5) = fx(1);
  pv1(k,6) = fy(1);
  
  xp(k,:) = x(:);
  yp(k,:) = y(:);
  
  # find the energies of the system
  pet = sum(pe);
  ket = sum(ke);
  
  time(k) = k*dt;
  energy(k,1)=pet;
  energy(k,2)=ket;
  energy(k,3) = pet+ket;
end

#
# plot the results
#

# positions
plot(xp(:,1),yp(:,1), ";Object 1;", xp(:,2), yp(:,2), ";Object 2;")
title("Positions of First 2 Objects")
xlabel ("X")
ylabel ("Y")

#plot(xp(:,:),yp(:,:))
#title("Positions of All Objects")
#xlabel ("X")
#ylabel ("Y")

# Energies
#plot(energy(:,1), ";Potential;", energy(:,2), ";Kinetic;", energy(:,3), ";Total;")
#xlabel("Iteration")
#ylabel("Engery")

%% Boris A Algorithm for Relativstic Charged Particle

%% Constants
tic
q = 1.0; % part charge
m = 1.0; % Electron mass 
c = 1.0; % Speed of light

%% Initial conditions

init_pos = [0.9, 0.0, 0.0];  % Initial position (at the +1/2 position) 
init_u   = [0.1, 0.0, 0.0];  % Initial velocity multiplied by gamma 

%% Time parameters
dt = pi / 1E2;      % Time step
total_time = 1.0E3; % total time to run simulation
num_steps = ceil(total_time/dt);  % Number of time steps

%% Preallocate arrays

pos      = zeros(num_steps, 3);
u        = zeros(num_steps, 3);
e_Field  = zeros(num_steps, 3);
b_Field  = zeros(num_steps, 3);

% MATLAB function "zeros(M,N)" create matrices of all zeros of size MxN

%% Initialize variables

pos(1,:) = init_pos; % position at the 1/2 steps
u(1,:)   = init_u;   % velocity times gamma at nn = 1

%% Iterative loop over "nn" time steps

for nn = 1:num_steps

%% find the current values for E and B fields at the n+1/2 value

pos_plusHalf = pos(nn,:); % a 1x3 vector for [x,y,z]
% find E field (nn+1/2)
eFactor = ( 0.01.* (pos_plusHalf(1).^2 + pos_plusHalf(2).^2 )^(-3/2));
e_Field(nn,:) = eFactor.*[ pos_plusHalf(1) , pos_plusHalf(2) , 0];
% find B field (nn+1/2)
b_Field(nn,:) = [ 0 , 0 , sqrt( pos_plusHalf(1).^2 + pos_plusHalf(2).^2 ) ];


%% 1st half step push of particle due to E field [Equation 3]
u_minus = u(nn,:) + (q/m)*(dt/2).*(e_Field(nn,:)); 

%% rotation piece - Boris (2 step rotation) [Equation 4] u_minus -> u_plus
% Boris A = [Eqs. 6, 7a, 8, 9]    

% compute magnitude of u_minus [needed for gamma_minus calculation]
u_minus_mag = ( u_minus(:,1)^2 + u_minus(:,2)^2 + u_minus(:,3)^2 )^(0.5);

% compute gamma at the u_minus 'velocity' [inline equation between eqs. 2 & 3]
gamma_minus = ( 1 + (u_minus_mag/c)^2 )^(0.5);

% magnitude of B field at n+1/2 
B_mag = sqrt( b_Field(nn,1).^2 + b_Field(nn,2).^2 + b_Field(nn,3).^2 );
% (used MATLAB function "sqrt()" to compute the square root)

% [Equation 6] phase angle 
theta = (q * dt * B_mag) / ( m * gamma_minus );


% [Equation 7a]
t_vector = ( tan(theta / 2) ).*( b_Field(nn,:) ./ B_mag );
 

% [Equation 8] step 1 of rotation u_minus -> u_prime 
u_prime = u_minus + cross( u_minus , t_vector ); 
% (used MATLAB function "cross(A,B)" to compute the cross product between A and B)

% needed for Equation 9
t_squared = dot( t_vector , t_vector );
% (used MATLAB function "dot(A,B)" to compute the dot product between A and B)

% [Equation 9] updating the final piece of the rotation u_prime -> u_plus 
u_plus = u_minus + ( 2 / ( 1 + t_squared ) ).*( cross( u_prime , t_vector ) );


  %% final E field half time step push [Equation 5]
u(nn+1,:) = u_plus + ( q * dt / ( 2 * m ) ) .* e_Field(nn,:);

%% Update the position [Equation 1]

% needed to compute gamma for equation 1
u_mag = norm( u(nn+1,:) );
% (used MATLAB function "norm(vector)" to compute the magnitude of a vector)

gamma = sqrt(1+(u_mag/c)^2);

% [Equation 1]
pos(nn+1,:) = dt.*u(nn+1,:)./gamma + pos(nn,:);

end % nn loop
toc

figure(1); 
plot3(pos(:,1),pos(:,2),pos(:,3))
xlabel('X Position');ylabel('Y Position');zlabel('Z Position');
title(strcat('Relativistic motion of electron using Boris A Method',', ',num2str(dt),'=dt'),Interpreter='latex')



clear all; clc; close all;
nu  = 0.25;
Vs  = 1500;                       % Shear wave velocity
rho = 2.3;                        % t/m^3 
loading_duration = 0.36;
ti = 0;
xi = -90-2.5;
tf = loading_duration;



% Node grid
element_size = 2.5;  
L1     = 90;
xwidth = L1 + element_size;
ywidth = element_size+ element_size;
depth  = 40;


% Node spacing
[X,Y,Z] = meshgrid(-xwidth:element_size:xwidth, -ywidth:element_size:ywidth, -depth:element_size:0);
coords = [X(:) Y(:), Z(:)];


% Define the cube boundaries
x_min = -L1 + eps;
x_max =  L1 - eps;
y_min = -element_size + eps;
y_max =  element_size - eps;
z_min = -depth + element_size  + eps;
z_max =  0      + eps;


% Create a logical index for points outside the specified cube
outside_cube = ~(coords(:, 1) > x_min & coords(:, 1) < x_max & ...
                 coords(:, 2) > y_min & coords(:, 2) < y_max & ...
                 coords(:, 3) > z_min & coords(:, 3) < z_max);


% Use logical indexing to keep only the points outside the cube
coords = coords(outside_cube, :);
amplitude    =      0.01;
type         =         11;
omega        =         [];
dt           =       1.25e-4;
N            =       2^12; 
t            = dt*[0:N-1];

%%
A_t  = amplitude * impulse_load(ti, loading_duration, t, type, omega); %[m]
disp = zeros(size(coords,1)*3, size(t,2));
acc  = zeros(size(coords,1)*3, size(t,2));
vel  = zeros(size(coords,1)*3, size(t,2));
for i=1:size(coords,1)
    XYnodes = coords(i,1:2:3);
    [u, v, ddu, ddv, Vr] = Rayleigh_wave(A_t, Vs, nu, rho, dt, XYnodes, t, xi);
    disp((i-1)*3+1,:) = u(:);
    disp((i-1)*3+3,:) = v(:);
    acc((i-1)*3+1,:)  = ddu(:);
    acc((i-1)*3+3,:)  = ddv(:);
end


% saving matricies
writematrix(disp,'disp.txt','Delimiter',',')
writematrix(acc ,'acc.txt','Delimiter',',')
writematrix(vel, 'vel.txt','Delimiter',',')
writematrix(coords,'coords.txt','Delimiter',',')

%%
% controling points
A_t  = amplitude * impulse_load(ti, loading_duration, t, type, omega);
XYnodes = [-xwidth 0; 0 0; xwidth 0;-xwidth -30; 0 -30; xwidth -30];
[u, v, ddu, ddv, Vr] = Rayleigh_wave(A_t, Vs, nu, rho, dt, XYnodes, t, xi);

writematrix([t;u],'disxpcontrol.txt','Delimiter',',')
writematrix([t;v],'disypcontrol.txt','Delimiter',',')
writematrix([t;ddu] ,'accxcontrol.txt','Delimiter',',')
writematrix([t;ddv] ,'accycontrol.txt','Delimiter',',')
% subplot(2,1,1)
% plot(t,u(1:1:3,:));
% legend("p1", "p2", "p3")
% ylabel("u [m]")
% y1=ylim 
% grid("on")
% xlim([0,0.5])
% subplot(2,1,2)
% plot(t,u(4:1:6,:));
% legend("p4", "p5", "p6")
% xlim([0,0.5])
% % ylim([-2e-3 2e-3]) 
% ylabel("u [m]")
% xlabel("Time [s]")
% grid("on")



function p = impulse_load(ti, tf, t, type, omega)
% IMPULSE_LOAD                        % Given the time vector and the
% definition of amplitude, initial and final time, this function returns the
% impulse load vector.
% 6.1.2020   (Month.Day.Year)       % Adriano Trono

% Definition of the load type...
% 1 : bell Po(1-cos(2*pi*t/to))
% 2 : step
% 3 : complete sine wave +- (single)
% 4 : positive half sine wave (single)
% 5 : descending triangle
% 6 : ascending triangle
% 7 : symmetric triangle
% 8 : positive step and negative step
% 9 : cosine up to pi/2
% 10: sine from ti
% 11: mexican hat

p = zeros(1, length(t));

% Variable associated with useful indices
idx_useful = find((t >= ti) & (t <= tf));

% Useful indices for the first half...
idx_useful1 = find((t >= ti) & (t <= (tf + ti) / 2));

% Useful indices for the second half...
idx_useful2 = find((t >= (tf + ti) / 2) & (t <= tf));

% Switch based on the load type
switch type
    case 1
        % Bell
        p(idx_useful) = 1 - cos(2 * pi * (t(idx_useful) - ti) / (tf - ti));
    case 2
        % Step
        p(idx_useful) = 1;
    case 3
        % Sine wave +-
        p(idx_useful) = sin(2 * pi * (t(idx_useful) - ti) / (tf - ti));
    case 4
        % Positive sine wave
        p(idx_useful) = sin(pi * (t(idx_useful) - ti) / (tf - ti));
    case 5
        % Descending triangle
        p(idx_useful) = 1 - (t(idx_useful) - ti) / (tf - ti);
    case 6
        % Ascending triangle
        p(idx_useful) = (t(idx_useful) - ti) / (tf - ti);
    case 7
        % Symmetric triangle
        p(idx_useful1) = 2 * (t(idx_useful1) - ti) / (tf - ti);
        p(idx_useful2) = 1 - 2 * (t(idx_useful2) - 0.5 * (ti + tf)) / (tf - ti);
    case 8
        % Positive step and negative step
        p(idx_useful1) = 1;
        p(idx_useful2) = -1;
    case 9
        % Cosine up to pi/2
        p(idx_useful) = cos(pi * (t(idx_useful) - ti) / (tf - ti) / 2);
    case 10
        % Sine from ti
        p(t > ti) = sin(omega * (t(t > ti) - ti));
    case 11
        % Mexican hat
        maxp = max(mexihat(-2 * pi, 2 * pi, length(idx_useful)));
        p(idx_useful) = mexihat(-2 * pi, 2 * pi, length(idx_useful))/maxp;
end

end



function [u, v, ddu, ddv, Vr] = Rayleigh_wave(A_t, Vs, nu, rho, dt, ...
    XYnodes, t, xi)
% Rayleigh_wave calculates the displacement vector in x and y directions due
% to the passage of a Rayleigh plane wave.
% A_t    : time signal to be considered,
% Vs     : shear wave velocity
% rho    : material density
% XYnodes: matrix with the first column as x positions and the second column as y positions
% (of the points for which displacements and velocities are to be calculated);
% xi     : location where the Rayleigh wave starts
% nu     : Poisson's ratio

% Duplicate the number of points to avoid capturing the wave restart...
N = 2 * length(t);

% Soil properties...
G   = Vs^2 * rho;                  % kPa
E   = 2 * (1 + nu) * G;
kk  = sqrt(2 * (1 - nu) / (1 - 2 * nu));  % relation between Vp and Vs
Vp  = kk * Vs; % P-wave velocity

% Roots are the roots of the polynomial in the auxiliary variable xi = (Vr/Vs)^2
% (taken from Graff's book page 325 book - 340 pdf)
rootss  = roots([1 -8 24-16*kk^-2 -16*(1-kk^-2)]);
Vr = real(rootss(find(rootss < 1)))^.5 * Vs;

% Frequency vector...
df           = 1 / (N) / dt;
f            = df * [0:N-1];

% Fourier transform of the amplitude...
A_f      = fft([A_t zeros(1, length(A_t))]);
U   = zeros(size(XYnodes,1), N);
V   = zeros(size(XYnodes,1), N);
ddU = zeros(size(XYnodes,1), N);
ddV = zeros(size(XYnodes,1), N);

for j = 2:N
    k = 2*pi*f(j)/Vr;
    s = sqrt(k^2 - (2*pi*f(j)/Vs)^2);
    q = sqrt(k^2 - (2*pi*f(j)/Vp)^2);
    
    for i = 1:size(XYnodes,1)
        z = abs(XYnodes(i,2));

        U(i,j) =  A_f(j)  * (-exp(-q*z) + ...
            2*q*s/k^2/(s^2/k^2+1)*exp(-s*z)) * exp(-complex(0,1)*k*(XYnodes(i,1) - xi))/(-1 + ...
            2*q*s/k^2/(s^2/k^2+1));
        
        
        V(i,j) =  A_f(j)* complex(0,1)*( 2*q/k/(s^2/k^2+1) * exp(-s*z) - q/k*exp(-q*z))...
            * exp(-complex(0,1)*k*(XYnodes(i,1)-xi))/(-1 + ...
            2*q*s/k^2/(s^2/k^2+1));

        
        ddU(i,j) = -(2*pi*f(j))^2*U(i,j);
        ddV(i,j) = -(2*pi*f(j))^2*V(i,j);
    end
end

u   = transpose(ifft(transpose(U), 'symmetric'));
u   = u(:, 1:N/2); 
v   = transpose(ifft(transpose(V), 'symmetric'));
v   = v(:, 1:N/2); 
ddu = transpose(ifft(transpose(ddU), 'symmetric'));
ddu = ddu(:, 1:N/2); 
ddv = transpose(ifft(transpose(ddV), 'symmetric'));
ddv = ddv(:, 1:N/2); 

end

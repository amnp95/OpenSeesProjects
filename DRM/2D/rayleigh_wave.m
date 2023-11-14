function [u, v, ddu, ddv, Vr] = rayleigh_wave(A_t, Vs, nu, rho, dt, ...
    XYnodes, t, xi)
% rayleigh_wave calculates the displacement vectors in the x and y directions
% resulting from the passage of a Rayleigh wave.
% A_t     : Time signal to be considered for wave propagation
% Vs      : Shear wave velocity
% rho     : Material density
% XYnodes : Matrix with the x-y positions of points where displacements and velocities are desired
% xi      : Starting point of the Rayleigh wave
% nu      : Poissons ratio

% Double the number of points to avoid wave start detection...
N = 2 * length(t);

% Soil properties...
G   = Vs^2 * rho;                  % kPa
E   = 2 * (1 + nu) * G;
kk  = sqrt(2 * (1 - nu) / (1 - 2 * nu));  % Relation between Vp and Vs
Vp  = kk * Vs; % Compressional wave velocity

% Calculate roots of the polynomial equation related to wave velocity ratio
% xi = (Vr/Vs)^2 (taken from Graffs book, page 325)
raices = roots([1 -8 24-16*kk^-2 -16*(1-kk^-2)]);
Vr = real(raices(find(raices < 1)))^.5 * Vs;

% Frequency vector...
df = 1 / (N) / dt;
f = df * [0:N-1];

% Fourier transform of the amplitude...
A_f = fft([A_t zeros(1, length(A_t))]);
U = zeros(length(XYnodes), N);
V = zeros(length(XYnodes), N);
ddU = zeros(length(XYnodes), N);
ddV = zeros(length(XYnodes), N);

for j = 2:N
    k = 2 * pi * f(j) / Vr;
    s = sqrt(k^2 - (2 * pi * f(j) / Vs)^2);
    q = sqrt(k^2 - (2 * pi * f(j) / Vp)^2);
    
    for i = 1:length(XYnodes)
        z = abs(XYnodes(i, 2));
        
        U(i, j) =  A_f(j) * (-exp(-q*z) + ...
            2 * q * s / k^2 / (s^2 / k^2 + 1) * exp(-s*z)) * exp(-1i * k * (XYnodes(i, 1) - xi)) / (-1 + ...
            2 * q * s / k^2 / (s^2 / k^2 + 1));
        
        V(i, j) =  A_f(j) * 1i * k * (2 * q / k / (s^2 / k^2 + 1) * exp(-s*z) - q / k * exp(-q*z))...
            * exp(-1i * k * (XYnodes(i, 1) - xi)) / (-1 + ...
            2 * q * s / k^2 / (s^2 / k^2 + 1));
        
        ddU(i, j) = -(2 * pi * f(j))^2 * U(i, j);
        ddV(i, j) = -(2 * pi * f(j))^2 * V(i, j);
    end
end

u = transpose(ifft(transpose(U), 'symmetric'));
u = u(:, 1:N/2); 
v = transpose(ifft(transpose(V), 'symmetric'));
v = v(:, 1:N/2); 
ddu = transpose(ifft(transpose(ddU), 'symmetric'));
ddu = ddu(:, 1:N/2); 
ddv = transpose(ifft(transpose(ddV), 'symmetric'));
ddv = ddv(:, 1:N/2); 
end

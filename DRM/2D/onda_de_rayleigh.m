function [u, v, ddu, ddv, Vr] = onda_de_rayleigh(A_t,Vs,nu,rho,dt,...
    XYnodes,t,xi)
% onda_de_rayleigh calcula el vector de desplazamientos en x y en y a partir
% del paso de una onda plana de Rayleigth
% A_t    : señal de tiempo que se considera que va a pasar,
% Vs     : velocidad de onda de corte
% rho    : densidad del material
% XYnodes: matriz con primera columna posición x segunda columna posición y
% (de los puntos cuyos desplazamientos y velocidades se quieren conocer);
% xi     : lugar donde comienza el paso de la onda de Rayleigh
% nu     : Módulo de Poisson

% se duplica la cantidad de puntos para que no se capte el reinicio de la
% onda...
N = 2*length(t);

% datos del suelo...
G   =   Vs^2*rho;                  % kPa
E   = 2*(1+nu)*G;
kk  = sqrt(2*(1-nu)/(1-2*nu));  %relacion entre Vp y Vs
Vp  = kk*Vs; % velocidad de onda P

% raices son las raices del polinomio en la variable auxiliar xi= (Vr/Vs)^2
%(sacado del libro de Graff pág 325 book - 340 pdf
raices  = roots([1 -8 24-16*kk^-2 -16*(1-kk^-2)]);
Vr = real(raices(find(raices<1)))^.5 * Vs;

% vector de frecuencias...
df           =     1/(N)/dt;
f            =   df*[0:N-1];

% transformada de Fourier de la amplitud...
A_f      = fft([A_t zeros(1,length(A_t))]);
U   = zeros(length(XYnodes),N);
V   = zeros(length(XYnodes),N);
ddU = zeros(length(XYnodes),N);
ddV = zeros(length(XYnodes),N);

for j = 2:N
    k = 2*pi*f(j)/Vr;
    s = sqrt(k^2 - (2*pi*f(j)/Vs)^2);
    q = sqrt(k^2 - (2*pi*f(j)/Vp)^2);
    
    for i = 1:length(XYnodes)
        z = abs(XYnodes(i,2));
        % se le agrega un menos adelante para que sea un movimiento
        % retrogradando con el avance
        % a las ecuaciones del Fede las multiplico por i
        % reocordá que V = - W en las filminas de Fede. 
%         U(i,j) = - A_f(j)* k  * (-exp(-q*z) + ...
%             2*q*s/k^2/(s^2/k^2+1)*exp(-s*z)) * exp(-complex(0,1)*k*(XYnodes(i,1) - xi));
%         V(i,j) = - A_f(j)* complex(0,1)*k*( 2*q/k/(s^2/k^2+1) * exp(-s*z) - q/k*exp(-q*z))...
%             * exp(-complex(0,1)*k*(XYnodes(i,1)-xi));

        U(i,j) =  A_f(j)  * (-exp(-q*z) + ...
            2*q*s/k^2/(s^2/k^2+1)*exp(-s*z)) * exp(-complex(0,1)*k*(XYnodes(i,1) - xi))/(-1 + ...
            2*q*s/k^2/(s^2/k^2+1));
        
        
        V(i,j) =  A_f(j)* complex(0,1)*( 2*q/k/(s^2/k^2+1) * exp(-s*z) - q/k*exp(-q*z))...
            * exp(-complex(0,1)*k*(XYnodes(i,1)-xi))/(-1 + ...
            2*q*s/k^2/(s^2/k^2+1));
%         U(i,j) = - A_f(j) *(-exp(-q*z) + ...
%             2*q*s/k^2/(s^2/k^2+1)*exp(-s*z)) * exp(-complex(0,1)*k*(XYnodes(i,1) - xi));
%         V(i,j) = - A_f(j) * complex(0,1)*( 2*q/k/(s^2/k^2+1) * exp(-s*z) - q/k*exp(-q*z))...
%             * exp(-complex(0,1)*k*(XYnodes(i,1)-xi));
        
        ddU(i,j) = -(2*pi*f(j))^2*U(i,j);
        ddV(i,j) = -(2*pi*f(j))^2*V(i,j);
    end
end

u   = transpose(ifft(transpose(U)  ,'symmetric'));
u   = u(:,1:N/2); 
v   = transpose(ifft(transpose(V)  ,'symmetric'));
v   = v(:,1:N/2); 
ddu = transpose(ifft(transpose(ddU),'symmetric'));
ddu = ddu(:,1:N/2); 
ddv = transpose(ifft(transpose(ddV),'symmetric'));
ddv = ddv(:,1:N/2); 

end


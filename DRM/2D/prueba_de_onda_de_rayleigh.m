
clear all
nu  = .35;
Vs  = 229;                       % velocidad de onda de corte
rho = 1.8;                       % t/m^3 

duracion_de_carga = 9e-2;
ti = 0;
xi = -10;
tf = duracion_de_carga;

% grilla de nudos..
element_size = .5;                          % paso entre nudos
x            = [-10:element_size:-9.5 9.5:element_size:10]; % abscisa
y            = -(0:element_size:10);                         % ordenada
aux_x        = ones(length(y),1)*x;
aux_y        = transpose(y)*ones(1,length(x));
XYnodes      = [aux_x(:) aux_y(:)];

x            = [-10:element_size:10]; % abscisa
y            = -(9.5:element_size:10);                         % ordenada
aux_x        = ones(length(y),1)*x;
aux_y        = transpose(y)*ones(1,length(x));

XYnodes      = union(XYnodes,[aux_x(:) aux_y(:)],'rows','stable');

amplitud     = 0.001;

tipo         =         11;
omega        =         [];
dt           =       1e-4;
N            =       2^12; 
t            = dt*[0:N-1];

A_t  = amplitud*carga_impulsiva(ti,duracion_de_carga,t,tipo,omega); %[m]

[u, v, ddu, ddv, Vr] = onda_de_rayleigh(A_t,Vs,nu,rho,dt,XYnodes,t,xi);




%%
scale = 500;
for i=1:6:length(t)-1
    clf 'reset'
    plot(XYnodes(:,1) + scale * u(:,i) , XYnodes(:,2) + scale * v(:,i),'.')
%         plot(XYnodes(:,1)  , XYnodes(:,2),'.')
    
        axis([-11 11 -11 1])
    shg
    hold on
    grid on
%     title(['t =' num2str(t(j)) '    carga = ' num2str(cargay(j))]);
    pause(dt/10000)
end
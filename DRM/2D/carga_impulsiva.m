function p = carga_impulsiva(ti,tf,t,tipo, omega)
% CARGA_IMPULSIVA                     % a partir del vector tiempo y de la
% definición de la amplitud y del tiempo inicial y final devuelve el vector
% de carga impulsiva...
% 6.1.2020   (Mes.Día.Año)       % Adriano Trono                    

% definición del tipo de carga...
% 1 : campana Po(1-cos(2*pi*t/to))
% 2 : escalón
% 3 : onda seno completa +- (una sola)
% 4 : media onda seno + (una sola)
% 5 : triángulo descendente
% 6 : triángulo creciente
% 7 : triángulo simétrico
% 8 : escalón positivo y escalón negativo
% 9 : coseno hasta pi/2
% 10: seno desde ti
% 11: mexican hat

p = zeros(1,length(t));

% variable asociada a los indices útiles 
idx_utiles = find((t>=ti)&(t<=tf));

% indices utiles primera mitad...
idx_utiles1 = find((t>=ti)&(t<=(tf+ti)/2));

% indices utiles segunda mitad...
idx_utiles2 = find((t>=(tf+ti)/2)&(t<=tf));

% indices utiles

switch tipo
    case 1
        % campana
        p(idx_utiles) = 1 - cos(2*pi*(t(idx_utiles)-ti)/(tf-ti));
    case 2
        % escalón
        p(idx_utiles) = 1;    
    case 3
        % onda seno +-
        p(idx_utiles) = sin(2*pi*(t(idx_utiles)-ti)/(tf-ti));
    case 4
        % onda seno +
        p(idx_utiles) = sin(pi*(t(idx_utiles)-ti)/(tf-ti));
    case 5
        % triángulo descendente
        p(idx_utiles) = 1 - (t(idx_utiles)-ti)/(tf-ti);
    case 6
        % triángulo creciente
        p(idx_utiles) = (t(idx_utiles)-ti)/(tf-ti);
    case 7
        % triángulo simétrico
        p(idx_utiles1) = 2 * (t(idx_utiles1)-ti)/(tf-ti);
        p(idx_utiles2) = 1 - 2*(t(idx_utiles2)-.5*(ti+tf))/(tf-ti);
    case 8
        % esalón positivo y escalón negativo
        p(idx_utiles1) =   1;
        p(idx_utiles2) = - 1;
    case 9 
        % coseno hasta pi/2
        p(idx_utiles) = cos(pi*(t(idx_utiles)-ti)/(tf-ti)/2);
    case 10
        % seno desde ti
        p(t>ti) = sin(omega*(t(t>ti)-ti));
    
    case 11
        % mexican hat
        p(idx_utiles) = mexihat(-2*pi, 2*pi, length(idx_utiles));
              
end
        
end






        
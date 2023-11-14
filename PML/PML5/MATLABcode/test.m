%testing functions
XYelement(1) = -13;
XYelement(2) = -12;
XYelement(3) = -13;
XYelement(4) = -13;
XYelement(5) = -12;
XYelement(6) = -13;
XYelement(7) = -12;
XYelement(8) = -12;

beta_0_x = 259.041;
beta_0_y = 259.041;
alpha_0_x = 1.61901;
alpha_0_y = 1.61901;
L_PML_x  = -8;
L_PML_y  = -8;
xi = -13;
yj = -13;
rho = 2000;
E = 2.08e+08;
nu = 0.3;
K = k_sym_discontinua5_PML22(XYelement,beta_0_x,beta_0_y,L_PML_x,L_PML_y,xi,yj,rho,E,nu);
M = m_sym_discontinua5_PML22(XYelement,alpha_0_x,alpha_0_y,L_PML_x,L_PML_y,xi,yj,rho,E,nu);
C = c_sym_discontinua5_PML2222(XYelement,alpha_0_x,alpha_0_y,beta_0_x,beta_0_y,L_PML_x,L_PML_y,xi,yj,rho,E,nu);

save('KMatlab.txt', 'K', '-ascii');
save('MMatlab.txt', 'M', '-ascii');
save('CMatlab.txt', 'C', '-ascii');




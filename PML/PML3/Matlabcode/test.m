%testing functions
XYelement(1) = 5;
XYelement(2) = 5;
XYelement(3) = 5;
XYelement(4) = 3;
XYelement(5) = 7;
XYelement(6) = 3;
XYelement(7) = 7;
XYelement(8) = 5;

beta_0_x = 2.0;
beta_0_y = 3.0;
alpha_0_x = 3.0;
alpha_0_y = 2.0;
L_PML_x  = 4;
L_PML_y  = 5;
xi = 5;
yj = 3;
rho = 2000;
E = 1.e4;
nu = 0.3;
K = k_sym_discontinua3_PML22(XYelement,beta_0_x,beta_0_y,L_PML_x,L_PML_y,xi,yj,rho,E,nu);
M = m_sym_discontinua3_PML22(XYelement,alpha_0_x,alpha_0_y,L_PML_x,L_PML_y,xi,yj,rho,E,nu);
C = c_sym_discontinua3_PML2222(XYelement,alpha_0_x,alpha_0_y,beta_0_x,beta_0_y,L_PML_x,L_PML_y,xi,yj,rho,E,nu);





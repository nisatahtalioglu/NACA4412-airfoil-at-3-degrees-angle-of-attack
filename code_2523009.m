clc;
clear;
alpha=3;
c=1;



C_p_u = readmatrix("naca4412@3-2.xlsx","Range" ,"C2 : C104");
C_p_l = readmatrix("naca4412@3-2.xlsx","Range", "C105 : C201");

x_u = readmatrix("naca4412@3-2.xlsx","Range", "A2 : A104");
x_l = readmatrix("naca4412@3-2.xlsx","Range" , "A105 : A201");

y_u = readmatrix("naca4412@3-2.xlsx","Range", "B2: B104");
y_l = readmatrix("naca4412@3-2.xlsx","Range", "B105: B201");

bound = transpose(linspace(0, 1, 96));    % integral boundaries
x_u = x_u(1:97,:);
x_l = x_l(1:97,:);
y_u = y_u(1:97,:);
y_l = y_l(1:97,:);

C_p_u = C_p_u(1:96,:);
C_p_l = C_p_l(1:96,:);


dy_u_dx_u = diff(y_u)./diff(x_u);
dy_l_dx_l = diff(y_l)./diff(x_l);

x_u = x_u(1:96,:);
x_l = x_l(1:96,:);
y_u = y_u(1:96,:);
y_l = y_l(1:96,:);

integral_1_1 = trapz(bound, C_p_l);
integral_1_2 = -trapz(bound, C_p_u);

integral_1 = integral_1_1 + integral_1_2;

C_n = (1/c) * integral_1;       % normal force coeff 
    
integral_2_1 = trapz(bound, C_p_u .* dy_u_dx_u);
integral_2_2 = trapz(bound, -C_p_l .* dy_l_dx_l);

integral_2 = integral_2_1 + integral_2_2;

C_a = (1/c) * integral_2;       % axial force coeff 
  
integral_3_1 = trapz(bound, C_p_l.*x_l);
integral_3_2 = trapz(bound,-C_p_u.*x_u);

integral_3 = integral_3_1 + integral_3_2;

integral_4_1 = trapz(bound,  C_p_u .* dy_u_dx_u.*y_u);
integral_4_2 = trapz(bound, -C_p_l .* dy_l_dx_l.*y_l);

integral_4 = integral_4_1 + integral_4_2;

C_m = 1/(c^2) * (integral_4 + integral_3);


C_d = C_n.*sind(alpha) + C_a.*cosd(alpha);

C_l = C_n.*cosd(alpha) - C_a.*sind(alpha);
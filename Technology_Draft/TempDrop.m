function output = TempDrop(h, T_lb_in, P, W, P_H_1, P_H_2, height, H_extr, LengthR, q2_0, c_p, r_Co, r_Ci, r_Go, r_Gi, lambda_F, lambda_c, lambda_G)

H = linspace(0, height, length(h));
z = linspace(-height/2, height/2, length(h));

T_sat = XSteam('Tsat_p', P);

%CLADDING
A = q2_0 * (P_H_1*LengthR + P_H_2*(1-LengthR)) * H_extr/(pi*W*c_p) * sin(pi*height/(2*H_extr)) + T_lb_in;
B = q2_0 * (P_H_1*LengthR + P_H_2*(1-LengthR)) * H_extr/(pi*W*c_p);

C_ci = q2_0 *((r_Co/lambda_c) * log(r_Co/r_Ci) + 1./h);
C_Fc = q2_0 *(r_Co/lambda_c *log(r_Co/r_Ci) + r_Co/lambda_G * log(r_Go/r_Gi) + r_Co/(2*lambda_F) + 1./h);

T_ci = A + B *sin(pi*z/H_extr) + C_ci.*cos(pi*z/H_extr);
T_Fc = A + B *sin(pi*z/H_extr) + C_Fc.*cos(pi*z/H_extr);

% T_Ci_max = A + sqrt(B^2 + C_ci.^2);
% z_Ci_max = H_extr/pi * atan(B./C_ci);
% 
% T_ci = temp(1, 1:length(temp));
% T_Fc = temp(2, 1:length(temp));

output = [T_ci;T_Fc];

end
function output = DFM_subc(T_lb, i_l, P_const, W, z_OSV, A_1, A_2, LengthR, height)
% SubCooled part

H = linspace(0, height, length(T_lb));
z = linspace(-height/2, height/2, length(T_lb));

T_sat = XSteam('Tsat_p', P_const);

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

g = 9.82; %m/s^2, earth gravity acc

rho_g = XSteam('rhoV_T', T_sat); %kg/m^3, saturated vapour phase density
rho_f = XSteam('rhoL_T', T_sat); %kg/m^3, saturated vapour phase density

G = W/(A_1*LengthR + A_2*(1-LengthR));

%Calculating the equilibrium quality for the whole length

for i = 1:length(H)
    x_e(i) = (i_l(i) - i_f)/i_fg;
end

% Surface tension calculations

for i = 1:length(z)
    sigma(i) = XSteam('st_T', T_lb(i)); %N/m, surface tension over the whole average channel
end

% Corrected calculation of x_a

for i = 1:length(H)
    x_a(i) = x_e(i) - x_e(find(H>=z_OSV,1))*exp(x_e(i) / x_e(find(H>=z_OSV,1)) - 1);
end

for i = 1:length(H)
    U_vj_s(i) = 2.9 * (sigma(i) * g * (rho_f - rho_g) / rho_f^2)^0.25;
end

% Corrected calculation of J_v_s and J_l_s
J_v_s = x_a * G / rho_g;
J_l_s = (1 - x_a) * G / rho_f;
J_s = J_v_s + J_l_s;

b = (rho_g / rho_f)^0.1;

for i = 1:length(H)
    beta(i) = 1 / (1 + (rho_g / rho_f) * (1 - x_a(i)) / x_a(i));
end

% Corrected calculation of C_0_s
for i = 1:length(H)
    C_0_s(i) = beta(i) * (1 + (1 / beta(i))^b);
end

% Corrected calculation of alpha_s
for i = 1:length(H)
    alpha_s(i) = J_v_s(i) / (C_0_s(i) * J_s(i) + U_vj_s(i));
end

output = alpha_s;

end
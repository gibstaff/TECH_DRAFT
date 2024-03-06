function output = Mult_r_2(T_lb, i_l, P_const, height)

H = linspace(0, height, length(i_l));
z = linspace(-height/2, height/2, length(i_l));

T_sat = XSteam('Tsat_p', P_const);

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

 %kg/m^3, saturated vapour phase density

for i = 1:length(T_lb)
    rho_g(i) = XSteam('rhoV_T', T_lb(i));
end

%kg/m^3, saturated liquid phase density

for i = 1:length(T_lb)
    rho_f(i) = XSteam('rhoL_T', T_lb(i));
end

% Equilibrium quality @ exit

x_e_ex = (i_l(length(i_l)) - i_f)/i_fg;

% Multiplier calculations

for i = 1:length(T_lb)
    r_2(i) = (rho_f(i)/rho_g(i)-1)*x_e_ex;
end

output = r_2;

end
function output = Mult_r_4(T_lb, i_l, alpha, P_const, height)

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

% Multiplier calculation

for i = 1:length(T_lb)
    r_4(i) = 1 - (rho_f(i) - rho_g(i))/rho_f(i) * alpha(i);
end

output = r_4;

end
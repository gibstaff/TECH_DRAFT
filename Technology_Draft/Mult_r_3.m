function output = Mult_r_3(T_lb, i_l, P_const, height)

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

%Dynamic Viscosity of the saturated fluid phase

mu_f = zeros(1,length(T_lb));
for i = 1:find(T_lb>=T_sat,1)-1
    mu_f(i) = XSteam('my_ph', P_const, XSteam('hL_T',T_lb(i))); %[N/m]
end
for i = find(T_lb>=T_sat,1):length(T_lb)
    mu_f(i) = XSteam('my_ph',P_const,XSteam('hL_P',P_const));
end

%[N/m] %Dynamic Viscosity of the saturated gas phase

mu_g = zeros(1,length(T_lb));
for i = 1:find(T_lb>=T_sat,1)-1
    mu_g(i) = XSteam('my_ph', P_const, XSteam('hV_T',T_lb(i))); %[N/m]
end
for i = find(T_lb>=T_sat,1):length(T_lb)
    mu_g(i) = XSteam('my_ph', P_const, XSteam('hV_P',P_const));
end

% Equilibrium quality

for i = 1:length(H)
    x_e(i) = (i_l(i) - i_f)/i_fg;
end

% Multiplier caluclation

for i = 1:length(x_e)
    r_3(i) = (1+(mu_f(i)/mu_g(i)-1)*x_e(i))^(-0.25) * (1+(rho_f(i)/rho_g(i)-1)*x_e(i));
end

output = r_3;

end
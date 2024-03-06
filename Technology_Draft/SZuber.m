function output = SZuber(T_lb, i_l, P_const, W, q2, D_h_1, D_h_2, A_1, A_2, LengthR, height)
% Using Saha-Zuber model to find the OSV point assuming constant pressure.



H = linspace(0, height, length(T_lb));
z = linspace(-height/2, height/2, length(T_lb));

T_sat = XSteam('Tsat_p', P_const);

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

%Re and Pr number calculations

for i = 1:find(T_lb >= T_sat, 1)-1
    Cp(i) = XSteam('Cp_pT', P_const, T_lb(i));
end

Cp = 1000*Cp; %J/(kg*Celsius) - liquid phase heat capacity in single-phase flow region

G = W/(A_1*LengthR + A_2*(1-LengthR)); %kg/(m^2*s), mass flow flux rate

for i = 1:find(T_lb >= T_sat, 1)-1
    mu(i) = XSteam('my_pT', P_const, T_lb(i)); %Pa*s, viscosity in the single-phase flow region
end

for i = 1:find(T_lb >= T_sat, 1)-1
    Re(i) = G*(D_h_1*LengthR + D_h_2*(1-LengthR))/mu(i);
end


% Calculating Prandtl number for the single-phase flow region

for i = 1:find(T_lb >= T_sat, 1)-1
    lambda(i) = XSteam('tc_pT', P_const, T_lb(i)); %W/(m*K), heat conductivity in the single-phase flow region
end

for i = 1:find(T_lb >= T_sat, 1)-1
    Pr(i) = Cp(i)*mu(i)/lambda(i);
end

% Saha-Zuber
%Let us define the Peclet number in the subcooled region

for i = 1:find(T_lb >= T_sat, 1)-1
    Pe(i) = Re(i)*Pr(i);
end

%Calculating the equilibrium quality for the whole length

for i = 1:length(H)
    x_e(i) = (i_l(i) - i_f)/i_fg;
end

for i = 1:length(Pe)
    x_OSV(i) = -154 * q2(i)/(G*i_fg);
end

z_OSV = H(find(x_e(1:length(x_OSV)) > x_OSV, 1));

output = [z_OSV; find(x_e(1:length(x_OSV)) > x_OSV, 1)];

end
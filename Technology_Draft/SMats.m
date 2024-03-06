function output = SMats(T_lb, P, D_h_1, D_h_2, A_1, A_2, LengthR, q2_avg, W, height)
% SMats - calculates the points of ONB and SUB using Sato-Matsumura
% correlation.
%
%T_lb is a vector of liquid bulk temperature, i_l is a vector of enthalpy
%in the coolant, P is the pressure at the inlet(assumed constant), D_h_1
%and D_h_2 are the hydraulic diameters, A_1 and A_2 are the flow areas,
%LengthR is the ratio of the partial rod length to the full one, q2_avg is
%a vector of the heat fluxes, W is the mass flow rate flux, hight is the
%active core height in m.


H = linspace(0, height, 1000);
z = linspace(-height/2, height/2, 1000);

T_sat = XSteam('Tsat_p', P);

for i = 1:find(T_lb == max(T_lb), 1)-1
    Cp(i) = XSteam('Cp_pT', P, T_lb(i));
end

Cp = 1000*Cp; %J/(kg*Celsius) - liquid phase heat capacity in single-phase flow region

G = W/(A_1*LengthR + A_2*(1-LengthR)); %kg/(m^2*s), mass flow flux rate

for i = 1:find(T_lb >= T_sat, 1)-1
    mu(i) = XSteam('my_pT', P, T_lb(i)); %Pa*s, viscosity in the single-phase flow region
end

for i = 1:find(T_lb >= T_sat, 1)-1
    Re(i) = G*(D_h_1*LengthR + D_h_2*(1-LengthR))/mu(i);
end


% Calculating Prandtl number for the single-phase flow region

for i = 1:find(T_lb >= T_sat, 1)-1
    lambda(i) = XSteam('tc_pT', P, T_lb(i)); %W/(m*K), heat conductivity in the single-phase flow region
end

for i = 1:find(T_lb >= T_sat, 1)-1
    Pr(i) = Cp(i)*mu(i)/lambda(i);
end


% Dittus and Boelter correlation is valid for Re>10e+4, 0.7<Pr<100,
% L/Dh>60. All is valid, so we use it here. n = 0.4, due to heating.

for i = 1:find(T_lb >= T_sat, 1)-1
    Nu(i) = 0.023 * Re(i)^0.8 * Pr(i)^0.4;
end

% As we use a rod bundle here, we use Markoczy correlation for a rod
% bundle.

B = 4/pi * ((1.295*10^-2)/(10.3*10^-3))^2 - 1;

for i = 1:find(T_lb >= T_sat, 1)-1
    Nu(i) = Nu(i) * (1+0.91*Re(i)^(-0.1)*Pr(i)^(0.4)*(1-2*exp(-B)));
end

%From the definition of Nusselt number, we find h.

for i = 1:find(T_lb >= T_sat, 1)-1
    h(i) = lambda(i)*Nu(i)/(D_h_1*LengthR + D_h_2*(1-LengthR));
end

% Sato-Matsumura, finding Tonb from Tsat + dTsuperheat

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

rho_g = XSteam('rhoV_T', T_sat); %kg/m^3, saturated vapour phase density

for i = 1:find(T_lb >= T_sat, 1)-1
    sigma(i) = XSteam('st_T', T_lb(i));
end

for i = 1:find(T_lb >= T_sat, 1)-1
    T_ONB(i) = T_sat + sqrt(8*sigma(i)*(T_sat+237.15)*q2_avg(i)/(lambda(i)*i_fg*rho_g));
end

for i = 1:find(T_lb >= T_sat, 1)-1
    T_w(i) = T_lb(i) + q2_avg(i)/h(i);
end

z_ONB = H(find(T_w >= T_ONB, 1)); %m, point of the ONB

z_SUB = H(find(T_lb == T_sat,1)); %m, point of the end of the subcooled boiling

output = [z_ONB, z_SUB; find(T_w > T_ONB, 1), find(T_lb == T_sat,1)];

end
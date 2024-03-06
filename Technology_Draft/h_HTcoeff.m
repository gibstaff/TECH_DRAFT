function output = h_HTcoeff(T_lb, q2, z_SUB, z_OSV, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, height)
% Calculates the convective heat transfer coefficient over the whole length of the
% assembly. P is in bar, q2 is in W/m^2. Output is in W/(m^2*K).

H = linspace(0, height, length(q2));
z = linspace(-height/2, height/2, length(q2));

T_sat = XSteam('Tsat_p', P);

h = zeros(1,length(T_lb));
T_w = zeros(1,length(T_lb));

%% Single-phase

for i = 1:find(H >= z_SUB, 1)-1
    Cp(i) = XSteam('Cp_ph', P, T_lb(i));
end

Cp = 1000*Cp; %J/(kg*Celsius) - liquid phase heat capacity in single-phase flow region

G = W/(A_1*LengthR + A_2*(1-LengthR)); %kg/(m^2*s), mass flow flux rate

for i = 1:find(H >= z_SUB, 1)-1
    mu(i) = XSteam('my_pT', P, T_lb(i)); %Pa*s, viscosity in the single-phase flow region
end

for i = 1:find(H >= z_SUB, 1)-1
    Re(i) = G*(D_h_1*LengthR + D_h_2*(1-LengthR))/mu(i);
end


% Calculating Prandtl number for the single-phase flow region

for i = 1:find(H >= z_SUB, 1)-1
    lambda(i) = XSteam('tc_pT', P, T_lb(i)); %W/(m*K), heat conductivity in the single-phase flow region
end

for i = 1:find(H >= z_SUB, 1)-1
    Pr(i) = Cp(i)*mu(i)/lambda(i);
end


% Dittus and Boelter correlation is valid for Re>10e+4, 0.7<Pr<100, L/Dh>60. All is valid, so we use it here. n = 0.4.

for i = 1:find(H >= z_SUB, 1)-1
    Nu(i) = 0.023 * Re(i)^0.8 * Pr(i)^0.4;
end

% As we use a rod bundle here, we use Markoczy correlation for a rod
% bundle.

B = 4/pi * ((1.295*10^-2)/(10.3*10^-3))^2 - 1;

for i = 1:find(H >= z_SUB, 1)-1
    Nu(i) = Nu(i) * (1+0.91*Re(i)^(-0.1)*Pr(i)^(0.4)*(1-2*exp(-B)));
end

% From the definition of Nusselt number, we find h and T_w.

for i = 1:find(H >= z_OSV, 1)-1
    h(i) = lambda(i)*Nu(i)/(D_h_1*LengthR + D_h_2*(1-LengthR));
    T_w(i) = T_lb(i)+q2(i)/h(i);
end

%% Two-phase Subcooled boiling

for i = 1:find(H >= z_SUB, 1)-1
    sigma(i) = XSteam('st_T', T_lb(i));
end

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

rho_g = XSteam('rhoV_T', T_sat); %kg/m^3, saturated vapour phase density


for i = find(H >= z_OSV, 1):find(T_lb >= T_sat, 1)-1
    % T_w(i) = 22.65*(q2(i)/10^6)^(0.5)*exp(-P/87) + T_sat; %Thom et al but
    % we have to use Sato-Matsumura here because we used it for ONB calcs
    T_w(i) = T_sat + sqrt(8*sigma(i)*(T_sat+237.15)*q2(i)/(lambda(i)*i_fg*rho_g));
    h(i) = q2(i)/(T_w(i)-T_lb(i));
end

%% Two-phase Saturated boiling

for i = find(T_lb >= T_sat, 1):length(T_lb)
    h(i) = 5.5 * (P/10)^(0.25) * (q2(i))^(2/3); %Rossohin estimation, seems very non-accurate
    T_w(i) = T_lb(i)+q2(i)/h(i);
end

output = [h;T_w];

end
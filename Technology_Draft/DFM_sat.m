
function output = DFM_sat(T_lb, i_l, P_const, W, A_1, A_2, D_h_1, D_h_2, LengthR, height)
%Gives the axial voidage distribution from saturated boiling using the DFM
%model. Output is a 1xlength(T_lb) vector.

H = linspace(0, height, length(T_lb));
z = linspace(-height/2, height/2, length(T_lb));

T_sat = XSteam('Tsat_p', P_const);

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

%Calculating the equilibrium quality for the whole length

for i = 1:length(H)
    x_e(i) = (i_l(i) - i_f)/i_fg;
end

G = W/(A_1 * LengthR + A_2 * (1-LengthR));

% DFM-model for saturated region

g = 9.82; %m/s^2, earth gravity acc

rho_g = XSteam('rhoV_T', T_sat); %kg/m^3, saturated vapour phase density
rho_f = XSteam('rhoL_T', T_sat); %kg/m^3, saturated vapour phase density

D_h = D_h_1 * LengthR + D_h_2 * (1-LengthR); %Averaged hydraulic diameter

for i = 1:length(z)
    sigma(i) = XSteam('st_T', T_lb(i)); %N/m, surface tension over the whole average channel
end

breakpoint = find(T_lb == T_sat,1);

for i = breakpoint:length(z)
    my = 9.13*10^-5; %obsobs
    J_l(i) = (1 - x_e(i)) * G / rho_f;
    J_v(i) = x_e(i) * G / rho_g;
    J(i) = J_v(i) + J_l(i);

    C_0_1 = 1.2;
    U_vj_1(i) = 1.41 * (sigma(i) * g * (rho_f - rho_g) / rho_f^2) ^(0.25);
    alpha_1(i) = J_v(i) / (C_0_1 * J(i) + U_vj_1(i));

    C_0_2 = 1.15;
    U_vj_2(i) = 0.35 * (D_h * g * (rho_f - rho_g) / rho_f) ^(0.5);
    alpha_2(i) = J_v(i) / (C_0_2 * J(i) + U_vj_2(i));

    C_0_3 = 1.05;
    U_vj_3(i) = 23 * (my * J_l(i) / (rho_g * D_h))^0.5 * (rho_f - rho_g) / rho_f;
    alpha_3(i) = J_v(i) / (C_0_3 * J(i) + U_vj_3(i));

    C_0_4 = 1;
    U_vj_4(i) = 1.53 * (sigma(i) * g * (rho_f - rho_g) / rho_g^2) ^(0.25);
    alpha_4(i) = J_v(i) / (C_0_4 * J(i) + U_vj_4(i));
end

alpha = zeros(size(alpha_1)); % Initialize alpha with zeros

for i = 1:length(alpha_1)
    if alpha_1(i) < 0.25 && alpha_1(i)>0
        alpha(i) = alpha_1(i);
    elseif alpha_2(i) < 0.75 && alpha_2(i) > 0.25
        alpha(i) = alpha_2(i);
    elseif alpha_3(i) < 0.95 && alpha_3(i) > 0.75
        alpha(i) = alpha_3(i);
    elseif alpha_4(i) < 1 && alpha_4(i) > 0.95
        alpha(i) = alpha_4(i);
    else
        alpha(i) = 0;
    end
end

output = alpha;

end
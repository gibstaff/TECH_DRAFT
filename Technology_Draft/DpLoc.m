function output = DpLoc(T_lb, q2, W, P_const, T_lb_in, A_1, A_2, P_H_1, P_H_2, LengthR, height, Spacer_num)

H = linspace(0, height, length(q2));
z = linspace(-height/2, height/2, length(q2));

T_sat = XSteam('Tsat_p', P_const);

i_g = XSteam('hV_T', T_sat)*10^3; %J/kg, saturated vapour phase enthalpy
i_f = XSteam('hL_T', T_sat)*10^3; %J/kg, saturated liquid phase enthalpy
i_fg = abs(i_f-i_g); %J/g, latent heat

i_l_in = XSteam('h_pT', P_const, T_lb_in); %kJ/kg, inlet enthalpy
i_l_in = i_l_in*1000; %J/kg

%kg/m^3, saturated vapor phase density

for i = 1:length(T_lb)
    rho_g(i) = XSteam('rhoV_T', T_lb(i));
end

%kg/m^3, saturated liquid phase density

for i = 1:length(T_lb)
    rho_f(i) = XSteam('rhoL_T', T_lb(i));
end

if length(W) ~= length(z)
    W = W*ones(1,length(z));
end

G = W./(A_1*LengthR + A_2*(1-LengthR));

% Spacers

d = height/Spacer_num; %Distance between spacers, assuming the final one is at the outlet [m]

sp1 = d;
sp2 = d + d;
sp3 = d + d*2;
sp4 = d + d*3;
sp5 = d + d*4;
sp6 = d + d*5;
sp7 = d + d*6;
sp8 = d + d*7;

% Local component

%1 note is 0,00318
sp_z1=sp1/0.00318;
sp_z2=sp2/0.00318;
sp_z3=sp3/0.00318;
sp_z4=sp4/0.00318;
sp_z5=sp5/0.00318;
sp_z6=sp6/0.00318;
sp_z7=sp7/0.00318;
sp_z8=sp8/0.00318;

i_sp1=i_l_in+(q2(150)*P_H_1)/W(150)*sp1;
i_sp2=i_l_in+(q2(300)*P_H_1)/W(300)*sp2;
i_sp3=i_l_in+(q2(450)*P_H_1)/W(450)*sp3;
i_sp4=i_l_in+(q2(600)*P_H_1)/W(600)*sp4;
i_sp5=i_l_in+(q2(749)*P_H_2)/W(749)*sp5;
i_sp6=i_l_in+(q2(899)*P_H_2)/W(899)*sp6;
i_sp7=i_l_in+(q2(1000)*P_H_2)/W(1000)*sp7;
i_sp8=i_l_in+(q2(1000)*P_H_2)/W(1000)*sp8;
%This inlet is wrong but idk if we want to do something about it

% Equilibrium quality calculations 

x_sp1 = (i_sp1 - i_f)/i_fg;
x_sp2 = (i_sp2 - i_f)/i_fg;
x_sp3 = (i_sp3 - i_f)/i_fg;
x_sp4 = (i_sp4 - i_f)/i_fg;
x_sp5 = (i_sp5 - i_f)/i_fg;
x_sp6 = (i_sp6 - i_f)/i_fg;
x_sp7 = (i_sp7 - i_f)/i_fg;
x_sp8 = (i_sp8 - i_f)/i_fg;

%Let us find two-phaser multipliers for each spacer using HEM

phi_lo1 = 1 + (rho_f(150)/rho_g(150) - 1)*x_sp1;
phi_lo2 = 1 + (rho_f(300)/rho_g(300) - 1)*x_sp2;
phi_lo3 = 1 + (rho_f(450)/rho_g(450) - 1)*x_sp3;
phi_lo4 = 1 + (rho_f(600)/rho_g(600) - 1)*x_sp4;
phi_lo5 = 1 + (rho_f(749)/rho_g(749) - 1)*x_sp5;
phi_lo6 = 1 + (rho_f(899)/rho_g(899) - 1)*x_sp6;
phi_lo7 = 1 + (rho_f(1000)/rho_g(1000) - 1)*x_sp7;
phi_lo8 = 1 + (rho_f(1000)/rho_g(1000) - 1)*x_sp8;

% The signle phase local loss coefficients are assumed to be

E_1 = 0.8;
E_2 = 0.8;


% The total local friction loss is then the sum of all the local losses
% and can be found as

for i = 1:length(z)
    if i> 0
        dp_local(i)= ((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo1);
    if i> 125*1
        dp_local(i)= dp_local(1)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo2);
    if i> 125*2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
        dp_local(i)= dp_local(126)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo3);
    if i> 125*3
        dp_local(i)= dp_local(125*2+1)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo4);
    if i> 125*4
        dp_local(i)= dp_local(125*3+1)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo5);
    if i> 125*5
        dp_local(i)= dp_local(125*4+1)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo6);
    if i> 125*6
        dp_local(i)= dp_local(125*5+1)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo7);
    if i> 127*7
        dp_local(i)= dp_local(125*6+1)+((G(i)^2)/(rho_f(i)*2))* E_1*(phi_lo8);
    end
    end
    end
    end
    end
    end
    end
    end

end

dp_local(1) = 0;
output = dp_local;

end
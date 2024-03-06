function output = Tlb_Enth(T_lb_in, W, P_const, LengthR, P_H_1, P_H_2, q2, height)
% Tlb_Enth - calculates the axial enthalpy and liquid bulk distribution for
% a channel with constant assumed pressure. The outputs are a vector of two
% rows, the channel is discretized into 1000 parts.
%
%T_lb_in - inlet liquid bulk temperature in Celsius, W - mass flow flux
%rate per assembly in kg/s, P_const - coolant pressure in bars, LengthR -
%ratio of partial to full length rod lengths, P_H_1 and P_H_2 - heated
%channel perimeters, q2 - axial heat flux distribution in W/m^2 (a vector),
%height - active core height.

z = linspace(-height/2, height/2, length(q2));

% Enthalpy destribution (coolant)

i_l_in = XSteam('h_pT', P_const, T_lb_in); %kJ/kg, inlet enthalpy
i_l_in = i_l_in*1000; %J/kg

% We take the weighted contribution of the heated perimeters

if length(W) ~= length(z)
    W = W*ones(1,length(z));
end

for i = 1:length(z)
    if i == 1
        i_l(i) = i_l_in + 1/W(i) * (q2(i) * (P_H_1*LengthR + P_H_2*(1-LengthR)) * z(i) - q2(i) * (P_H_1*LengthR + P_H_2*(1-LengthR)) * z(i));
    end
    if i > 1
        i_l(i) = i_l(i-1) + 1/W(i) * (q2(i) * (P_H_1*LengthR + P_H_2*(1-LengthR)) * z(i) - q2(i) * (P_H_1*LengthR + P_H_2*(1-LengthR)) * z(i-1));
    end
end


% Liquid bulk temp calculation

for i = 1:length(i_l)
    T_lb(i) = XSteam('T_ph',P_const,i_l(i)/1000); %Liquid bulk temp distr
end

output = [T_lb; i_l];

end
function output = DpFric(T_lb, r_3, P_const, W, A_1, A_2, D_h_1, D_h_2, P_w_1, P_w_2, LengthR, dRod, dWRod, height)

H = linspace(0, height, length(T_lb));
z = linspace(-height/2, height/2, length(T_lb));

T_sat = XSteam('Tsat_p', P_const); %Celsius

%Dynamic Viscosity of the saturated fluid phase

mu_f = zeros(1,length(T_lb));
for i = 1:find(T_lb>=T_sat,1)-1
    mu_f(i) = XSteam('my_ph', P_const, XSteam('hL_T',T_lb(i))); %[N/m]
end
for i = find(T_lb>=T_sat,1)-1:length(T_lb)
    mu_f(i) = XSteam('my_ph',P_const,XSteam('hL_P',P_const));
end

% Mass flow flux

G_1 = W./A_1;
G_2 = W./A_2;

if length(W) ~= length(z)
    W = W*ones(1,length(T_lb));
end

G = W./(A_1*LengthR + A_2*(1-LengthR));
D_h = D_h_1*LengthR + D_h_2*(1-LengthR);

% Reynolds number

Re = zeros(1,length(T_lb));
for i = 1:find(H>=2.134,1)-1
    Re(i) = (G_1*D_h_1)/mu_f(i);
end
for i = find(H>=2.134,1):length(T_lb)
    Re(i) = (G_2*D_h_2)/mu_f(i);
end

% for i = 1:length(T_lb)
%     Re(i) = (G(i)*D_h)/mu_f(i);
% end

%kg/m^3, saturated liquid phase density

for i = 1:length(T_lb)
    rho_f(i) = XSteam('rhoL_T', T_lb(i));
end

k=0.003; %Assumption of friction coefficient

% Fanning coefficient calculation

Cf = zeros(1,length(T_lb));
% for i = 1:find(H>=2.134,1)-1
%     Cf(i) = (-3.6*log10(((k/D_h_1)/3.7)^1.11 + 6.9/Re(i)))^-2;
% end
% for i = find(H>=2.134,1):length(T_lb)
%     Cf(i) = (-3.6*log10(((k/D_h_2)/3.7)^1.11 + 6.9/Re(i)))^-2;
% end

%Alyoshin's correlation which shouldn't be applicable here according to the
%source (but is fine according to the compendium lmao)
for i = 1:find(H>=2.134,1)-1
    Cf(i) = 0.38 * (P_w_1/(92*2*pi*dRod/2 + 2*2*pi*dWRod/2)) * (A_1/(92*pi*(dRod/2)^2 + 2*pi*(dWRod/2)^2))^(0.45) * Re(i)^(-0.25);
    %
end
for i = find(H>=2.134,1):length(T_lb)
    Cf(i) = 0.38 * (P_w_2/((92-14)*2*pi*dRod/2 + 2*2*pi*dWRod/2)) * (A_1/((92-14)*pi*(dRod/2)^2 + 2*pi*(dWRod/2)^2))^(0.45) * Re(i)^(-0.25);
    %
end

% for i = 1:length(T_lb)
%     Cf(i) = (-3.6*log10(((k/D_h)/3.7)^1.11 + 6.9/Re(i)))^-2;
% end

%Calculating frictional pressure loss

dp_fric = zeros(1, length(T_lb));
for i = 1:find(H>=2.134,1)-1
    dp_fric(i) = r_3(i)*(G_1^2/(2*rho_f(i))) * height * (4*Cf(i)/D_h_1); %[Pa]
end
for i = find(H>=2.134,1):length(T_lb)
    dp_fric(i) = r_3(i)*(G_2^2/(2*rho_f(i))) * height * (4*Cf(i)/D_h_2); %[Pa]
end

% for i = 1:length(T_lb)
%     dp_fric(i) = r_3(i)*(G(i)^2/(2*rho_f(i))) * height * (4*Cf(i)/D_h);
% end

dp_fric(1) = 0;

output = dp_fric;

end
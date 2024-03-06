clear

CorePower = 3926; %MWth, the utilized thermal power
CorePower = CorePower*10^6; %W
FuelPower = CorePower*0.962; %W, the heat deposited into the fuel rods+cladding - percentage taken from THydraulics slides

CoreEqD = 5.163; %m - equivalent diameter of a "circle" made of equal and equispaced squares (asemblies)
CoreH = 3.81; %m - core active height

CoreV = 3.81 * pi * (5.163/2)^2; %m^3

AvgFuelP = FuelPower/CoreV; %W/m^3 %Core-mean power density


% Geometric calculations, contd

w = 13.80; %cm, channel lattice dimension
w = 13.80/100; %m

dRod = 10.3; %mm, outer rod diameter
dRod = dRod*10^-3; %m, outer rod diameter

dWRod = 25.21; %mm, outer water rod diameter
dWRod = dWRod*10^-3; %m, outer water rod diameter

P_H_1 = 92*pi*dRod; %m, heated perimeter for both full and part length rod region
P_H_2 = (92-14)*pi*dRod; %m, heated perimeter for full length rod region

P_w_1 = 4*w + 92*pi*dRod + 2*pi*dWRod;
P_w_2 = 4*w + (92-14)*pi*dRod + 2*pi*dWRod;

D_H_1 = (4*w^2 - 92*pi*dRod^2 - 2*pi*dWRod^2)/P_H_1; %m, heated diameter for both full and part length rod region
D_H_2 = (4*w^2 - (92-14)*pi*dRod^2 - 2*pi*dWRod^2)/P_H_2; %m, heated diameter for full length rod region

A_1 = P_H_1*D_H_1/4; %m, flow area for both full and part length rod region
A_2 = P_H_2*D_H_2/4; %m, heated diameter for full length rod region

D_h_1 = (4*w^2 - 92*pi*dRod^2 - 2*pi*dWRod^2)/P_w_1; %m
D_h_2 = (4*w^2 - (92-14)*pi*dRod^2 - 2*pi*dWRod^2)/P_w_2; %m

A_rods = P_H_1*2.134 + P_H_2*(3.81-2.134); %m^2, outer surface A of fuel rods in one assembly

% Axial peaking factor calculation

%Equation is from THydraulics exercise 01_P02

R_extr = 6*(CoreEqD/2)/5; %m, extrapolation radius (reflection), R/Rextr = 5/6 = H/Hextr
H_extr = 6*CoreH/5; %m, extrapolation height

fP_axial = (pi*CoreH)/(2*H_extr*sin((pi*CoreH)/(2*H_extr)));
fP_radial = (2.405*CoreEqD/2)/(2*R_extr*besselj(1, 2.405*(CoreEqD/2)/R_extr));


% Bundle-average power calculation

%Let's start by calculating average power per volume

q3_coreAvg = CorePower/CoreV; %W/m^3, average power density per core volume

V_ass = 3.81*w^2; %m^3, the total active volume of an assembly

q_faAvg = fP_axial*q3_coreAvg*V_ass; %W, average power of an average assembly

while abs(q_faAvg*872 - CorePower) > 10e-6
    q_faAvg = q_faAvg*CorePower/(q_faAvg*872);
end

q3_coreAvg = q_faAvg/fP_radial/V_ass;
q_faHot = fP_radial*fP_axial*q3_coreAvg*V_ass; %W, average power of a hot assembly

q2_0_avg = q_faAvg/A_rods; %W/m^2
q2_0_hot = q_faHot/A_rods; %W/m^2

LengthR = 2.134/3.810;


% Modelling heat flux distribution

z = linspace(-CoreH/2,CoreH/2,1000);
H = linspace(0,CoreH,1000);

q2_avg = q2_0_avg*cos(pi*z/H_extr);
q2_hot = q2_0_hot*cos(pi*z/H_extr);


% Enthalpy destribution (coolant) and liquid bulk distribution

T_lb_in = 278; %Celsius

W = 14502; %kg/s, value from IAEA report
W = W/872; %kg/s, recalculated for flow for subassembly

P = 7.07; %MPa
P = P*10; %bar - the pressure will need to be iterated to account for the drop to improve the model. For now, let's assume constant.

temp = Tlb_Enth(T_lb_in, W, P, LengthR, P_H_1, P_H_2, q2_avg, CoreH);

T_lb = temp(1,1:length(q2_avg));
i_l = temp(2, 1:length(q2_avg));

T_sat = XSteam('Tsat_p', P);
dT_subc_in = T_sat-T_lb_in;


% Sato-Matsumura, finding Tonb from Tsat + dTsuperheat

clear temp

temp = SMats(T_lb, P, D_h_1, D_h_2, A_1, A_2, LengthR, q2_avg, W, CoreH);

z_ONB = temp(1,1)
z_SUB = temp(1,2)

% Saha-Zuber, finding the significant void start

clear temp
temp = SZuber(T_lb, i_l, P, W, q2_avg, D_h_1, D_h_2, A_1, A_2, LengthR, CoreH);

z_OSV = temp(1,1)

% Hot channel calculations for enthalpy and temperature

clear temp
temp = Tlb_Enth(T_lb_in, W, P, LengthR, P_H_1, P_H_2, q2_hot, CoreH);

T_lb_h = temp(1, 1:length(q2_hot));
i_l_h = temp(2, 1:length(q2_hot));

% Hot channel Sato-Matsumura correlation to find zONB and zSUB

clear temp
temp = SMats(T_lb_h, P, D_h_1, D_h_2, A_1, A_2, LengthR, q2_hot, W, CoreH);

z_ONB_h = temp(1,1)
z_SUB_h = temp(1,2)

% Hot channel Saha-Zuber model to find zOSV

clear temp
temp = SZuber(T_lb_h, i_l_h, P, W, q2_hot, D_h_1, D_h_2, A_1, A_2, LengthR, CoreH);

z_OSV_h = temp(1,1)


% DFM-model for saturated region, average channel

alpha = DFM_sat(T_lb, i_l, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);

% DFM-model for subcooled region, average channel

alpha_s = DFM_subc(T_lb, i_l, P, W, z_OSV, A_1, A_2, LengthR, CoreH);

% Might be prudent to continue Saha-Zuber until it is exceeded by saturated
% boiling to evade the drop at the saturation region start? If not, remove
% this part.

subset_indices = find(H>=z_OSV,1):find(H>=z_SUB,1) + find(alpha(find(H>=z_SUB,1):length(H)) >= alpha_s(find(H>=z_SUB,1):length(H)),1);

alpha(subset_indices) = alpha_s(subset_indices);

% DFM-model for sat region, hot channel

alpha_h = DFM_sat(T_lb_h, i_l_h, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);

% DFM-model for subc region, hot channel

alpha_s_h = DFM_subc(T_lb_h, i_l_h, P, W, z_OSV_h, A_1, A_2, LengthR, CoreH);

subset_indices2 = find(H>=z_OSV_h,1):find(H>=z_SUB_h,1) + find(alpha_h(find(H>=z_SUB_h,1):length(H)) >= alpha_s_h(find(H>=z_SUB_h,1):length(H)),1);

alpha_h(subset_indices2) = alpha_s_h(subset_indices2);
% Pressure drop, initial data

% Multipliers

% r_2=19;
% r_3 = 8;
% r_3 = Mult_r_3(i_l,P,CoreH);
% r_4=0.15;

G = linspace(W/(A_1*LengthR + A_2*(1-LengthR))*0.5, W/(A_1*LengthR + A_2*(1-LengthR))*1.5, 1000);

% rho_f = XSteam('rhoL_T', T_sat); %kg/m^3, saturated vapour phase density


% DpFric(T_lb, r_3, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);
% DpLoc(q2_avg,W,P,T_lb_in,A_1,A_2,P_H_1,P_H_2,LengthR,CoreH,8);
%% Average channel flow characteristic

% When multiplying the average assembly power by some constant, the heat
% flux for that power is essentially the nominal power multiplied by that
% constant. So, we can just use the multiplied flux in calculations.

W2 = G.*(A_1*LengthR + A_2*(1-LengthR)); %kg/s

% dp_grav = r_4*H(length(H))*rho_f; %Doesn't depend on G or W.

for i = 2:6
    clear temp Tlbtemp ilbtemp
    temp = Tlb_Enth(T_lb_in, W2, P, LengthR, P_H_1, P_H_2, q2_avg*(0.25*i),CoreH);
    Tlbtemp = temp(1,1:length(temp));
    ilbtemp = temp(2,1:length(temp));

    clear r_3 r_2
    r_3 = Mult_r_3(Tlbtemp, ilbtemp, P, CoreH);
    r_2 = Mult_r_2(Tlbtemp, ilbtemp, P, CoreH);
    
    clear rho_f
    for l = 1:length(temp(1,1:length(temp)))
        rho_f(l) = XSteam('rhoL_T', temp(1,l));
    end

    for j = 1:length(W2)
        
        clear temp1 temp_OSV temp_s
        temp1 = DFM_sat(Tlbtemp, ilbtemp, P, W2(j), A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);
        temp_OSV = SZuber(Tlbtemp, ilbtemp, P, W2(j), q2_avg.*(0.25*i), D_h_1, D_h_2, A_1, A_2, LengthR, CoreH);
        temp_OSV = temp_OSV(1);
        temp_s = DFM_subc(Tlbtemp, ilbtemp, P, W2(j), temp_OSV, A_1, A_2, LengthR, CoreH);

        temp_indice = find(H>=temp_OSV(1,1),1):find(temp(1:length(temp))>=T_sat,1) + find(alpha(find(temp(1:length(temp))>=T_sat,1):length(H)) >= alpha_s(find(temp(1:length(temp))>=T_sat,1):length(H)),1);

        temp_1(temp_indice) = temp_s(temp_indice);

        r_4 = Mult_r_4(temp(1,1:length(temp)),temp(2,1:length(temp)),temp1,P,CoreH);

        clear temp2
        temp2 = DpFric(temp(1,1:length(temp)), r_3, P, W2(j), A_1, A_2, D_h_1, D_h_2, P_w_1, P_w_2, LengthR, dRod, dWRod, CoreH);
        dp_fric(i-1,j) = temp2(length(temp2));

        clear temp2
        temp2 = DpLoc(temp(1,1:length(temp)),q2_avg*(0.25*i),W2(j),P,T_lb_in,A_1,A_2,P_H_1,P_H_2,LengthR,CoreH,8);
        dp_local(i-1,j) = temp2(length(temp2));

        clear temp2
        temp2 = r_2.*(G(j)^2)./rho_f;
        dp_acc(i-1,j) = temp2(length(temp2));

        clear temp2
        temp2 = r_4.*CoreH.*rho_f;
        dp_grav(i-1,j) = temp2(length(temp2));

        dp(i-1,j) = dp_fric(i-1,j) + dp_local(i-1,j) + dp_acc(i-1,j) + dp_grav(i-1,j);
    end
end

%% Plots

find(G>=1.526717629825628e+03,1)

figure
plot(G,-dp(1,1:length(dp))/10^5)
hold on
plot(G,-dp(2,1:length(dp))/10^5)
plot(G,-dp(3,1:length(dp))/10^5)
plot(G,-dp(4,1:length(dp))/10^5)
plot(G,-dp(5,1:length(dp))/10^5)
hold off
xlim([760, 2300])
title('Outlet Pressure vs. G for different Q')

% Adjust the size of the figure
fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
fontsize(gca, 8, "points")

xlabel('G [$\frac{kg}{m^{2}s}$]', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('$\Delta$P [bar]', 'Interpreter', 'latex', 'FontSize', 10)

legend({'0.5 x Q','0.75 x Q','1 x Q','1.25 x Q','1.5 x Q'},'Location','southwest', 'FontSize', 6)

saveas(gcf, 'Pressure_G.pdf');
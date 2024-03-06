clear

CorePower = 3926; %MWth, the utilized thermal power
CorePower = CorePower*10^6; %W
FuelPower = CorePower*0.962; %W, the heat deposited into the fuel rods+cladding - percentage taken from THydraulics slides

CoreEqD = 5.163; %m - equivalent diameter of a "circle" made of equal and equispaced squares (asemblies)
CoreH = 3.81; %m - core active height

CoreV = 3.81 * pi * (5.163/2)^2; %m^3

AvgFuelP = FuelPower/CoreV; %W/m^3 %Core-mean power density


%% Geometric calculations, contd

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

%% Axial peaking factor calculation

%Equation is from THydraulics exercise 01_P02

R_extr = 6*(CoreEqD/2)/5; %m, extrapolation radius (reflection), R/Rextr = 5/6 = H/Hextr
H_extr = 6*CoreH/5; %m, extrapolation height

fP_axial = (pi*CoreH)/(2*H_extr*sin((pi*CoreH)/(2*H_extr)));
fP_radial = (2.405*CoreEqD/2)/(2*R_extr*besselj(1, 2.405*(CoreEqD/2)/R_extr));


%% Bundle-average power calculation

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


%% Modelling heat flux distribution

z = linspace(-CoreH/2,CoreH/2,1000);
H = linspace(0,CoreH,1000);

q2_avg = q2_0_avg*cos(pi*z/H_extr);
q2_hot = q2_0_hot*cos(pi*z/H_extr);


%% Enthalpy destribution (coolant) and liquid bulk distribution

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


%% Sato-Matsumura, finding Tonb from Tsat + dTsuperheat

clear temp

temp = SMats(T_lb, P, D_h_1, D_h_2, A_1, A_2, LengthR, q2_avg, W, CoreH);

z_ONB = temp(1,1)
z_SUB = temp(1,2)

%% Saha-Zuber, finding the significant void start

clear temp
temp = SZuber(T_lb, i_l, P, W, q2_avg, D_h_1, D_h_2, A_1, A_2, LengthR, CoreH);

z_OSV = temp(1,1)

%% Hot channel calculations for enthalpy and temperature

clear temp
temp = Tlb_Enth(T_lb_in, W, P, LengthR, P_H_1, P_H_2, q2_hot, CoreH);

T_lb_h = temp(1, 1:length(q2_hot));
i_l_h = temp(2, 1:length(q2_hot));

%% Hot channel Sato-Matsumura correlation to find zONB and zSUB

clear temp
temp = SMats(T_lb_h, P, D_h_1, D_h_2, A_1, A_2, LengthR, q2_hot, W, CoreH);

z_ONB_h = temp(1,1)
z_SUB_h = temp(1,2)

%% Hot channel Saha-Zuber model to find zOSV

clear temp
temp = SZuber(T_lb_h, i_l_h, P, W, q2_hot, D_h_1, D_h_2, A_1, A_2, LengthR, CoreH);

z_OSV_h = temp(1,1)


%% DFM-model for saturated region, average channel

alpha = DFM_sat(T_lb, i_l, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);

%% DFM-model for subcooled region, average channel

alpha_s = DFM_subc(T_lb, i_l, P, W, z_OSV, A_1, A_2, LengthR, CoreH);

% Might be prudent to continue Saha-Zuber until it is exceeded by saturated
% boiling to evade the drop at the saturation region start? If not, remove
% this part.

subset_indices = find(H>=z_OSV,1):find(H>=z_SUB,1) + find(alpha(find(H>=z_SUB,1):length(H)) >= alpha_s(find(H>=z_SUB,1):length(H)),1);

alpha(subset_indices) = alpha_s(subset_indices);

%% DFM-model for sat region, hot channel

alpha_h = DFM_sat(T_lb_h, i_l_h, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);

%% DFM-model for subc region, hot channel

alpha_s_h = DFM_subc(T_lb_h, i_l_h, P, W, z_OSV_h, A_1, A_2, LengthR, CoreH);

subset_indices2 = find(H>=z_OSV_h,1):find(H>=z_SUB_h,1) + find(alpha_h(find(H>=z_SUB_h,1):length(H)) >= alpha_s_h(find(H>=z_SUB_h,1):length(H)),1);

alpha_h(subset_indices2) = alpha_s_h(subset_indices2);

%% Pressure drop, initial data

% Multipliers

% r_2=19;
% r_3 = 8;
% r_4 = 0.15;
r_2 = Mult_r_2(T_lb, i_l, P, CoreH);
r_3 = Mult_r_3(T_lb, i_l, P, CoreH);
r_4 = Mult_r_4(T_lb, i_l, alpha, P, CoreH);

G = W/(A_1*LengthR + A_2*(1-LengthR));

%kg/m^3, saturated vapour phase density

% rho_f = XSteam('rhoL_T', T_sat);
for i = 1:length(T_lb)
    rho_f(i) = XSteam('rhoL_T', T_lb(i));
end


%% Friction Component

dp_fric = DpFric(T_lb, r_3, P, W, A_1, A_2, D_h_1, D_h_2, P_w_1, P_w_2, LengthR, dRod, dWRod, CoreH);


%% Gravity component

%The value of sin(phi) is 1 for a vertical channel, since sin(90)=1.
%dp_grav = rho_f*(z_SUB*g+r_4*(L-z_SUB))); %[Pa]

dp_grav = r_4.*CoreH.*rho_f; %[Pa]


%% Acceleration component

dp_acc=r_2.*(G^2)./rho_f; %[Pa]


%% Local component

dp_local = DpLoc(T_lb, q2_avg, W, P, T_lb_in, A_1, A_2, P_H_1, P_H_2, LengthR, CoreH, 8);


%% Inlet and outlet pressure drop

zeta_in = 0.5;
dp_inlet = zeta_in*G*G/(2*rho_f(1));

zeta_rodchange = (1-A_1/A_2)^2;

dp_rodchange = zeta_rodchange*(W/A_1)^2/(2*rho_f(find(H>=2.134,1)));

zeta_ex = 1;

dp_exit = zeta_ex*G^2/(2*rho_f(length(rho_f))); %[Pa]?

%% Total drop, average channel

% Finally, we find the total pressure drop as the sum of all the
% components

dp = -dp_fric - dp_grav - dp_acc - dp_local;  %dp_rodchange; %[Pa]

dp(1) = dp(1) - dp_inlet;

dp(length(dp)) = dp(length(dp)) - dp_exit;

dp(find(H>=2.134,1)) = dp(find(H>=2.134,1)) - dp_rodchange;

for i = 1:length(dp)
    if i == 1
        P_test(i) = P;
    end
    if i == 2
        P_test(i) = P + dp(1)/10^5;
    end
    if i > 2
        P_test(i) = P_test(2) + dp(i)/10^5;
    end
end
 
%% HOT CHANNEL 

% r_3_h = r_3;
r_2_h = Mult_r_2(T_lb, i_l, P, CoreH);
r_3_h = Mult_r_3(T_lb, i_l, P, CoreH);
r_4_h = Mult_r_4(T_lb, i_l, alpha_h, P, CoreH);

dp_fric_h = DpFric(T_lb_h, r_3_h, P, W, A_1, A_2, D_h_1, D_h_2, P_w_1, P_w_2, LengthR, dRod, dWRod, CoreH);

dp_local_h = DpLoc(T_lb_h, q2_hot, W, P, T_lb_in, A_1, A_2, P_H_1, P_H_2, LengthR, CoreH, 8);

dp_acc_h = r_2_h.*(G^2)./rho_f; %[Pa]

dp_grav_h = r_4_h.*CoreH.*rho_f;

dp_h = -dp_fric_h - dp_grav_h - dp_acc_h - dp_local_h;  %dp_rodchange; %[Pa]

dp_h(1) = dp_h(1) - dp_inlet;

dp_h(length(dp_h)) = dp_h(length(dp_h)) - dp_exit;

dp_h(find(H>=2.134,1)) = dp_h(find(H>=2.134,1)) - dp_rodchange;

for i = 1:length(dp_h)
    if i == 1
        P_test_h(i) = P;
    end
    if i == 2
        P_test_h(i) = P + dp_h(1)/10^5;
    end
    if i > 2
        P_test_h(i) = P_test_h(2) + dp_h(i)/10^5;
    end
end


%% Convective heat transfer coefficient

clear temp
temp = h_HTcoeff(T_lb_h, q2_hot, z_SUB_h, z_OSV_h, P, W, A_1, A_2, D_h_1, D_h_2, LengthR, CoreH);

h = temp(1,1:length(temp));
T_w = temp(2,1:length(temp));


%% Fuel temp distr

%say aluminium
c_p = 5.414592721 *10^3;
lambda_c = 0.18*10^2; % cladding thermal conductivity
lambda_G = 0.3496 *10^3; %gas gap thermal conductivity;
lambda_F = (0.042+2.71*10^-4 *670)^-1 + 6.9*10^-11 *670^3;

r_Gi = 8.76*10^-3 /2;
r_Go = 4.82*10^-3;
r_Ci = r_Go;
r_Co = 10.3*10^-3;

clear temp
temp = TempDrop(h, T_lb_in, P, W, P_H_1, P_H_2, CoreH, H_extr, LengthR, q2_0_hot, c_p, r_Co, r_Ci, r_Go, r_Gi, lambda_F, lambda_c, lambda_G);

T_ci = temp(1, 1:length(temp));
T_Fc = temp(2, 1:length(temp));


%% Plots

% Heat flux

figure
plot(H, q2_avg)
hold on
plot(H, q2_hot)
hold off

fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
fontsize(gca, 8, "points")

xlabel('z [m$]$', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('q$^{"}$[$\frac{W}{m^2}$]', 'Interpreter', 'latex', 'FontSize', 10)

legend({'Average assembly', 'Hot assembly'}, 'Location', 'northeast', 'FontSize', 6)

saveas(gcf, 'HFlux_comp.pdf');

% Pressure

figure
subplot(2,1,1)
plot(H, P_test)
hold on
plot(H, P_test_h)
hold off

fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
fontsize(gca, 8, "points")

xlabel('z [m$]$', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('P$_{tot}$[bar]', 'Interpreter', 'latex', 'FontSize', 10)

legend({'Average assembly', 'Hot assembly'}, 'Location', 'northeast', 'FontSize', 6)
xlim([-0.1, 3.9])

subplot(2,1,2)
plot(H, (dp-dp_h)/10^5)

fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
fontsize(gca, 8, "points")

xlabel('z [m$]$', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('P$_{tot}$(hot) - P$_{tot}$(avg) [bar]', 'Interpreter', 'latex', 'FontSize', 10)
xlim([-0.1, 3.9])
 
saveas(gcf, 'Ptot_comp.pdf');


% Pressure 2

figure
plot(H, P_test)
hold on
plot(H, P_test_h)
hold off

fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
fontsize(gca, 8, "points")

xlabel('z [m$]$', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('P$_{tot}$[bar]', 'Interpreter', 'latex', 'FontSize', 10)

legend({'Average assembly', 'Hot assembly'}, 'Location', 'northeast', 'FontSize', 6)
xlim([-0.1, 3.9])

saveas(gcf, 'Ptot_comp2.pdf');


% Voidage

alpha(subset_indices) = alpha_s(subset_indices);
alpha_h(subset_indices2) = alpha_s_h(subset_indices2);

figure 
plot(H, alpha)
hold on
% plot(H(subset_indices), alpha_s(subset_indices))
plot(H, alpha_h)
% plot(H(subset_indices2), alpha_s_h(subset_indices2))
hold off

%Adjust the size of the figure
fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
fontsize(gca, 8, "points")

xlabel('z [m$]$', 'Interpreter', 'latex', 'FontSize', 10)
ylabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 10)

legend({'Average assembly', 'Hot assembly'}, 'Location', 'southeast', 'FontSize', 6)

saveas(gcf, 'Alpha_comp.pdf');


figure
plot(z,T_ci)
xlabel('z [m]');
ylabel('Temperature [C]');
%legend('T(max) = 589C','Location','northwest')
%legend('z(max) = 0.14m','Location','northwest')

%Adjust the size of the figure
fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
saveas(gcf, 'Cladding_dist.pdf');


figure
plot(z,T_Fc)

xlabel('z [m]');
ylabel('Temperature [C]');
%Adjust the size of the figure
fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
saveas(gcf, 'Fuel_dist.pdf');


figure
plot(z,h)

xlabel('z [m]');
ylabel('h');
%Adjust the size of the figure
fig = gcf; % Get the current figure handle
fig.PaperPositionMode = 'auto'; % Use the paper's position for the figure
fig_pos = fig.PaperPosition; % Get the current figure's position
fig.PaperSize = [fig_pos(3) fig_pos(4)]; % Set the paper size to match the figure size
saveas(gcf, 'HT_coeff.pdf');


%% Test


% % Testing saturation temp change with pressure
% 
% for i = 1:length(P_test)
%     test(i) = XSteam('Tsat_p', P_test_h(i));
% end
% 
% disp(test(1) - test(length(test)))
% 
% % Recalculating wall distrib.
% 
% i_l2 = q2_0_avg*(P_H_1*LengthR + P_H_2*(1-LengthR)) * H_extr/(W*pi) * (sin(pi*z/H_extr) + sin(pi*CoreH/(2*H_extr))) + XSteam('h_pT', P, T_lb_in)*1000;
% 
% for i = 1:find(H>=z_SUB_h)-1
%     Cp(i) = XSteam('Cp_ph', P, i_l_h(i)/1000);
% end
% 
% for i = find(H>=z_SUB_h):length(i_l_h)
%     Cp(i) = XSteam('CpL_T', T_sat);
% end
% 
% for i = 1:length(z)
%     T_lb2(i) = q2_0_hot*(P_H_1*LengthR + P_H_2*(1-LengthR))*H_extr/(W*pi*1000*Cp(i)) * (sin(pi*z(i)/H_extr) + sin(pi*CoreH/(2*H_extr))) + T_lb_in;
% end
% 
% for i = 1:length(z)
%     T_w2(i) = q2_0_hot*(P_H_1*LengthR + P_H_2*(1-LengthR))*H_extr/(W*pi*1000*Cp(i)) * (sin(pi*z(i)/H_extr) + sin(pi*CoreH/(2*H_extr))) + q2_0_hot/h(i)*cos(pi*z(i)/H_extr) + T_lb_in;
% end
% 
% for i = 1:length(z)
%     T_ci2(i) = q2_0_hot*(P_H_1*LengthR + P_H_2*(1-LengthR))*H_extr/(W*pi*1000*Cp(i)) * (sin(pi*z(i)/H_extr) + sin(pi*CoreH/(2*H_extr))) + q2_0_hot*(r_Co/lambda_c * log(r_Co/r_Ci) + 1/h(i))*cos(pi*z(i)/H_extr) + T_lb_in;
% end
% 
% plot(H, T_lb2)
% hold on
% plot(H, T_lb_h)
% plot(H, T_w)
% plot(H, T_w2)
% % plot(H, T_ci2)
% hold off
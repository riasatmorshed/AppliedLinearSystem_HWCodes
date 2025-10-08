clear; clc; close all;
A = [0 1;-26 -10];
B = [0 0;468 52];
% Typically, for a SISO system, the output y=x1, which corresponds to the
% first state variable, not the identity matrix for C
C = [1 0];
D = []; %since it's a SISO system
ss_obj = ss(A,B,C,D);
transferFunc = tf(ss_obj);
% (b)
k = dcgain(ss_obj);
% (c)
% DC gain is the magnitude of the frequency response at omega = 0
DC_gain = freqresp(ss_obj,0);
% (d)
[wn,zeta,p] = damp(ss_obj);
% (e)

polesFromTF = roots(double(cell2sym(transferFunc.Denominator(1))));
zerosFromTF = roots(double(cell2sym(transferFunc.Numerator(1))));


wn = abs(polesFromTF);  % Natural frequencies (rad/s)
zeta = -real(polesFromTF) ./ wn;  % Damping ratios

% Convert natural frequencies to Hz (from rad/s)
wn_Hz = wn / (2 * pi);

% Compute time constants
timeConstants = -1 ./ real(polesFromTF);  % Time constants

% Display the results
polesFromTF, zerosFromTF, wn_Hz, zeta, timeConstants
%% Problem 2 (My solve)

freq_rangeHz = [abs(polesFromTF);abs(zerosFromTF)] ./ (2*pi);
minFreqHz_lim = min(freq_rangeHz)/10;
maxFreqHz_lim = max(freq_rangeHz)*100;

freqHz = logspace(log10(minFreqHz_lim),log10(maxFreqHz_lim),200);
H = squeeze(freqresp(transferFunc,2*pi*freqHz));
dB_all = 20.* log10(H);
deg_all = rad2deg(angle(H));

poles_Hz = abs(polesFromTF) / (2 * pi);
zeros_Hz = abs(zerosFromTF) / (2 * pi);

figure()
subplot(2,1,1);
semilogx(freqHz, dB_all, 'LineWidth', 2);
hold on;
semilogx(freqHz, ones(1, length(freqHz)) * -3, 'k--', 'LineWidth', 2);  % Corrected the plot of -3 dB line
ylabel('dB Magnitude', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Log Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
title('Bode Plot for Input 1 (u)', 'FontSize',10,'FontWeight','bold')
xline(poles_Hz, 'r--', 'LineWidth', 1.5);  % Corrected xline syntax for poles
grid on;

subplot(2,1,2);
semilogx(freqHz, deg_all, 'LineWidth', 2);
ylabel('Phase (degrees)', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Log Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
title('Bode Plot for Input 2 ($\dot{u}$)', 'FontSize',10,'FontWeight','bold')
xline(poles_Hz, 'r--', 'LineWidth', 1.5);  % Corrected xline syntax for poles
grid on;

%% steady state comparison
U = 1;
jwHz = 8;
time = linspace(1,2.5,1000);
H_jw = squeeze(freqresp(transferFunc,2*pi*jwHz));
yp_1 = U * abs(H_jw(1)) * sin(2*pi*jwHz * time + angle(H_jw(1)));
yp_2 = U * abs(H_jw(2)) * sin(2*pi*jwHz * time + angle(H_jw(2)));
% by lsim function
[y_lsim_1,time] = lsim(ss_obj(1,1),sin(2*pi*jwHz * time),time); % I have two inputs so i need to make it single output
[y_lsim_2,time] = lsim(ss_obj(1,2),sin(2*pi*jwHz * time),time);
figure()
plot(time(:,1), yp_1(1,:), 'LineWidth',2);
hold on
plot(time(:,1),y_lsim_1(:,1), 'LineWidth',2)
ylabel('Response Comparison', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
title('Output for 1st Input (u)','FontSize', 10, 'FontWeight', 'bold')
legend('Harmonic Response','LSIM Output');
figure()
plot(time(:,1),yp_2(1,:), 'LineWidth',2);
hold on 
plot(time(:,1),y_lsim_2(:,1), 'LineWidth',2);
grid on
ylabel('Response Comparison', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
title('Output for 2nd Input ($\dot{u}$)','FontSize', 10, 'FontWeight', 'bold')
legend('Harmonic Response','LSIM Output');
grid on
%% Problem 03

% settling time
% Given values
zeta = 0.8;
omega_n = 6.25;
s1 = -5 + 6.25*sqrt(1-0.8^2);
s2 = -5 - 6.25*sqrt(1-0.8^2);
char_poly = poly([s1, s2]);

[wn, zeta_damped, poles_out] = damp(char_poly);


%% Problem 05

% (a)
A_5 = [-4 -1; 26 -6];
B_5 = [1; 0];
C_5 = [];
D_5 = [];

ss_obj_5 = ss(A_5,B_5,C_5,D_5);
[wn_5, zeta_damped_5, poles_out_5] = damp(ss_obj_5);

% (b)

s1_5 = -wn_5(1)*2 + 1j *wn_5(1)*2 * sqrt(1-0.8839^2)
s2_5 = -wn_5(1)*2 - 1j *wn_5(1)*2 * sqrt(1-0.8839^2)

syms s g1 g2 
G = [g1 g2];
sI_A_BG = s*eye(2) - ss_obj_5.a + ss_obj_5.b*G;
CE_5 = collect(det(sI_A_BG),s);
CED_5 = poly([s1_5,s2_5]);
LHS_5 = [1 0; 6 26];
RHS_5 = [CED_5(2)-10;CED_5(3)-50];
G_5 = (inv(LHS_5)*RHS_5)';
% checking the result
eig(ss_obj_5.a-ss_obj_5.b*G_5)
%% Problem 06

A_6 = [-10 18 0;-9 8 0;-9 5 3];
B_6 = [1;1;0];
C_6 = [];
D_6 = [];

ss_obj_6 = ss(A_6,B_6,C_6,D_6);
[wn_6, zeta_damped_6, poles_out_6] = damp(ss_obj_6);
 
wn_6_cl = 20;
zeta_6_cl = 7.2443 * zeta_damped_6(2);

s1_6 = -wn_6_cl + 1j * wn_6_cl * sqrt(1-zeta_6_cl^2)
s2_6 = -wn_6_cl - 1j * wn_6_cl * sqrt(1-zeta_6_cl^2)
s3_6 = -wn_6_cl*zeta_6_cl

syms s g1 g2 g3
G_6 = [g1 g2 g3];
sI_A_BG_6 = s*eye(3) - ss_obj_6.a + ss_obj_6.b*G_6;
CE_6 = collect(det(sI_A_BG_6),s);
CED_6 = poly([s1_6,s2_6,s3_6])
LHS_6 = [1 1 0; 7 -2 -4; -30 -3 -85]; %hardcoded
RHS_6 = [CED_6(2)+1;CED_6(3)-76;CED_6(4)+246];
G_6 = (inv(LHS_6)*RHS_6)'
% checking the result
eig(ss_obj_6.a-ss_obj_6.b*G_6)
%% Problem 07

A_7 = [0 1 0;0 0 1;-52 -30 -4];
B_7 = [0;0;1];
C_7 = [20 1 0];
D_7 = 0;

ss_obj_7 = ss(A_7,B_7,C_7,D_7);
[wn_7, zeta_damped_7, poles_out_7] = damp(ss_obj_7);

s1_7 = -200;
s2_7 = -600;
s3_7 = -1000;

syms s g1 g2 g3
G_7 = [g1 g2 g3];
sI_A_BG_7 = s*eye(3) - ss_obj_7.a + ss_obj_7.b*G_7;
CE_7 = collect(det(sI_A_BG_7),s)
CED_7 = poly([s1_7,s2_7,s3_7])
LHS_7 = [0 0 1;0 1 0;1 0 0]; %hardcoded
RHS_7 = [CED_7(2)-4;CED_7(3)-30;CED_7(4)-52];
G_7 = (inv(LHS_7)*RHS_7)'
% checking the result
eig(ss_obj_7.a-ss_obj_7.b*G_7)
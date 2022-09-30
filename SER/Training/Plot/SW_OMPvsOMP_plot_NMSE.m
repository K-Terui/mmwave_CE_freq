%% graph drawing phase
%% load data
clear
load("Theory.mat");
load("result_SW_OMP_M80.mat");
SW_OMP_NMSE = NMSE_dB;
load("result_OMP_M80.mat");
OMP_NMSE = NMSE_dB;
clear NMSE_dB;
load("result_SW_OMP_meanPL_M80.mat");
meanPL_SW_OMP_NMSE = NMSE_dB_meanPL;
load("result_OMP_meanPL_M80.mat");
meanPL_OMP_NMSE = NMSE_dB_meanPL;
clear NMSE_dB_meanPL;

%% SW-OMP vs OMP without theoretical value
figure
grid on
hold on

plot(sn_dB, OMP_NMSE, "r-s", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, SW_OMP_NMSE, "b-o", "LineWidth", 2, "MarkerSize", 10)


xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["OMP", "SW-OMP"], "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])

%% SW-OMP vs OMP with theoretical value
figure
grid on
hold on

plot(sn_dB, OMP_NMSE, "r-s", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, NMSE_theory_OMP, "r--*", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, SW_OMP_NMSE, "b-o", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, NMSE_theory, "b--*", "LineWidth", 2, "MarkerSize", 10)

xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["OMP(Simu)", "OMP(Conv)", "SW-OMP(Simu)", "SW-OMP(Conv)"], "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])

%% SW-OMP vs OMP without theoretical value (meanPL)
figure
grid on
hold on

plot(sn_dB, meanPL_OMP_NMSE, "r-s", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, meanPL_SW_OMP_NMSE, "b-o", "LineWidth", 2, "MarkerSize", 10)


xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["OMP", "SW-OMP"], "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])

%% SW-OMP vs OMP with theoretical value (meanPL)
figure
grid on
hold on

plot(sn_dB, meanPL_OMP_NMSE, "r-s", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, NMSE_theory_OMP, "r--s", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, meanPL_SW_OMP_NMSE, "b-o", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, NMSE_theory, "b--o", "LineWidth", 2, "MarkerSize", 10)


xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["OMP(Simu)", "OMP(Conv)", "SW-OMP(Simu)", "SW-OMP(Conv)"], "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])
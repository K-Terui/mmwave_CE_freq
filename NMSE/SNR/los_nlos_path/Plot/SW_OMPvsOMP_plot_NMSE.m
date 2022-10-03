%% graph drawing phase
clear
load("result_SW_OMP_LosNLos.mat");
figure('Position', [100 100 500 375]);
grid on
hold on

plot(sn_dB, NMSE_dB_losnlos, "r-s", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, NMSE_dB, "b-o", "LineWidth", 2, "MarkerSize", 10)

xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["w/ distinction","w/o distinction"], "Location","northeast")
% title(['NMSE versus SNR (M = ', num2str(M), ')'])
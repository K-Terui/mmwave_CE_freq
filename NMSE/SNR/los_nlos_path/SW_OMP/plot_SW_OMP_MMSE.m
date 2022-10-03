%% plot a glaph 
sn_dB = -15:5:10;         % SNR roop iterate (-15,-10,-5,0,5,10)[dB]
M = 100;
load("Theory.mat");
load("SW_OMP_NMSEdB");
load("MMSE_NMSEdB");
figure
grid on
hold on

plot(sn_dB, SW_OMP_NMSEdB, "b-o", "LineWidth", 2, "MarkerSize", 10)
plot(sn_dB, MMSE_NMSE_dB, "g-s", "LineWidth", 2, "MarkerSize", 10)

xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["SW-OMP", "MMSE-NMSE"], "Location","northeast")
title(['NMSE versus SNR (M = ', num2str(M), ')'])
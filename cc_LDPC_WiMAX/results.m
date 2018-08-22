% varying Nit,  fixed n = 576 and r = 1
clear all;
close all;

load results/n576N10.mat;
Pbit1 = Pbit;
SNR_dB1 = SNR_dB;

load results/n576N100.mat;
Pbit2 = Pbit;
SNR_dB2 = SNR_dB;

load results/n1344N10.mat;
Pbit3 = Pbit;
SNR_dB3 = SNR_dB;

load results/n1344N100.mat;
Pbit4 = Pbit;
SNR_dB4 = SNR_dB;


%Uncoded BER
SNR_dBU = [-1 : 2 , 2.5];
SNR = 10.^(SNR_dBU/10);
Pbit_uncoded = qfunc(sqrt(2*SNR));

% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNR_dB1,Pbit1,'k-',SNR_dB2,Pbit2,'rs-',SNR_dB3,Pbit3,'go-',SNR_dB4,Pbit4,'b+-',SNR_dBU,Pbit_uncoded,'b--','LineWidth',2,'MarkerSize',10)
axis([-1 2.5 1e-3 1e0])
hleg = legend('n = 576 ; N = 10',...
              'n = 576 ; N = 100',...
              'n = 1344 ; N = 10',...
              'n = 1344 ; N = 100','Uncoded BER');
set(hleg,'position',[0.15 0.13 0.20 0.15]);
xlabel('SNR $\Gamma$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
        'YGrid', 'on', 'XGrid', 'on');
    
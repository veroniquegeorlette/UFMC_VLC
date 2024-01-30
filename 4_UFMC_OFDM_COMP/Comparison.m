clear all;
close all;

OFDM_juxt_BER=load('OFDM_juxt.mat', 'BER');
OFDM_juxt_SNR=load('OFDM_juxt.mat', 'SNRplot');
UFMC_BER=load('myMat_DC.mat', 'BER1');
UFMC_SNR=load('myMat_DC.mat', 'SNRplot');


figure;
semilogy(OFDM_juxt_SNR.SNRplot,OFDM_juxt_BER.BER,'b-*','LineWidth',2,'DisplayName','OFDM');
hold on;
semilogy(UFMC_SNR.SNRplot,UFMC_BER.BER1,'r-o','LineWidth',2,'DisplayName','UFMC');
%hold on;
%semilogy(SNR_sans_DC.SNRplot,BER_sans_DC.BER1,'k-<','LineWidth',2,'DisplayName','Real without DC');
grid on;
xlabel('Signal to Noise Ratio');
ylabel('Bit Error Rate');
title ('BPSK OFDM/UFMC for real DCO');
legend;
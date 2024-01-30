clear all;
close all;
clc;

ufmc_ber_=load('without_bias.mat','BER1');
ufmc_snr_=load('without_bias.mat','SNRplot');
ufmc_dc_ber=load('with05.mat','BER1');
ufmc_dc_snr=load('with05.mat','SNRplot');


figure;

semilogy(ufmc_dc_snr.SNRplot,ufmc_dc_ber.BER1,'b-*','LineWidth',2,'DisplayName','With Bias');
hold on;
semilogy(ufmc_snr_.SNRplot,ufmc_ber_.BER1,'r-o','LineWidth',2,'DisplayName','Without DC Bias');

title ('Bit Error Rate vs SNR for UFMC without DC bias and with DC bias of 0.5');
grid on;
legend;
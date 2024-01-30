clear all;
close all;
clc;

OFDM_HS_BE=load('OFDM_norm_BER.mat','be');
OFDM_HS_SNR=load('OFDM_norm_SNR.mat','sn');
OFDM_juxt_BER=load('OFDM_DC_bias_BER.mat','be');
OFDM_juxt_SNR=load('OFDM_DC_bias_SNR.mat','sn');
OFDM_QPSK_BER=load('OFDM_05_BER.mat', 'be');
OFDM_QPSK_SNR=load('OFDM_05_SNR.mat', 'sn');
%OFDM_16QAM_BER=load('OFDM_16QAM_juxt.mat', 'BER');
%OFDM_16QAM_SNR=load('OFDM_16QAM_juxt.mat', 'SNRplot');


figure;

semilogy(OFDM_juxt_SNR.sn,OFDM_juxt_BER.be,'b-*','LineWidth',2,'DisplayName','Without Bias');
hold on;
semilogy(OFDM_HS_SNR.sn,OFDM_HS_BE.be,'r-o','LineWidth',2,'DisplayName','With DC Bias');
hold on;
semilogy(OFDM_QPSK_SNR.sn,OFDM_QPSK_BER.be,'k-<','LineWidth',2,'DisplayName','Bias of 0.5');
%grid on;
title ('Bit Error Rate vs SNR for OFDM without DC bias and with DC bias');
grid on;
legend;


clear all;
close all;
clc;

OFDM_BPSK_BER=load('UFMC_BPSK_juxt.mat', 'BER1');
OFDM_BPSK_SNR=load('UFMC_BPSK_juxt.mat', 'SNRplot');
OFDM_QPSK_BER=load('UFMC_QPSK_juxt.mat', 'BER1');
OFDM_QPSK_SNR=load('UFMC_QPSK_juxt.mat', 'SNRplot');
OFDM_16QAM_BER=load('UFMC_16QAM_juxt.mat', 'BER1');
OFDM_16QAM_SNR=load('UFMC_16QAM_juxt.mat', 'SNRplot');


figure;
semilogy(OFDM_BPSK_SNR.SNRplot,OFDM_BPSK_BER.BER1,'b-*','LineWidth',2,'DisplayName','BPSK');
hold on;
semilogy(OFDM_QPSK_SNR.SNRplot,OFDM_QPSK_BER.BER1,'r-o','LineWidth',2,'DisplayName','QPSK');
hold on;
semilogy(OFDM_16QAM_SNR.SNRplot,OFDM_16QAM_BER.BER1,'k-<','LineWidth',2,'DisplayName','16QAM');
grid on;
title ('Bit Error Rate vs SNR BPSK, QPSK and 16QAM for UFMC');
legend;


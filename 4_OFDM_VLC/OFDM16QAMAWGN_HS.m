tic
clear all;
close all;
bps=1;
Weather_impact=0.8;
nsymconst=2^bps;
NFFT=256;
GI=1/32;
nbsymbol=1000;%10000;
data= randi ([0 1],NFFT*nbsymbol*bps,1);
%vect=[8 9 13 12 10 11 15 14 2 3 7 6 0 1 5 4];
Qammodulator=comm.RectangularQAMModulator('ModulationOrder', nsymconst, ...
    'BitInput',true,...
    'NormalizationMethod', 'Average power'); 
      %  'SymbolMapping','Custom',...
    %'CustomSymbolMapping',vect,...
dataoutmapping=step(Qammodulator,data);

%Adding Hermitian symmetry to see if it is the same

Separate_symb_HS=reshape(dataoutmapping,nbsymbol,NFFT);

%data_HS_total=zeros(2*NFFT+2,nbsymbol);
for i=1:nbsymbol
    data_HS_total(i,:)=[0, Separate_symb_HS(i,:),0,fliplr(conj(Separate_symb_HS(i,:)))]; 
    dataOFDM(i,:)=ifft(data_HS_total(i,:),2*NFFT+2);
end

%%
dataOFDM_=dataOFDM';
dataifft=dataOFDM_([(end-(NFFT*GI-1)):end 1:end],:);
Bias=0.5;
DC_OFDM_symbols=dataifft+Bias;
DC_OFDM_symb_weather=Weather_impact*DC_OFDM_symbols;
%%
%Conversion parallelle à serie de tous les symboles
dataifft1=DC_OFDM_symb_weather(:);
% !!!!!! enlever le prefixe cyclique puis faire l'fft pour avoir en domaine
% freq
BER=[];
SNRplot=[];
for SNR=0:1:30
%calcul de la puissance moyenne du signal reçu
p = sum(dataifft1.*conj(dataifft1))/length(dataifft1);
snr = 10^(SNR/10);%valeur linéaire du rapport signal à bruit
desv = sqrt(p/snr)/sqrt(2);
n = randn(size(dataifft1));% %+ 1i*randn(size(dataifft1));
n = n*desv;
signalnoisy= dataifft1+ n; 
 
%converti en vecteur colonne, NB colonne=NB symboles
signalnoisyy=reshape(signalnoisy,(2*NFFT+2)+NFFT*GI,[]);
%Retire valeur DC
zeroDCOFDMrec=signalnoisyy-Bias;
%retire le Guard interval
zerosOFDMrealsR1=signalnoisyy(NFFT*GI+1:(2*NFFT+2)+NFFT*GI,:);
%Passage en domaine frequentiel
datafft=fft(zerosOFDMrealsR1,(2*NFFT+2)); 
%on ne prend que la premiere moitiee de la sortie de la FFT (HS)
Temp=datafft(2:NFFT+1,:)';
dataaffft=Temp(:);

variance=10^(-SNR/10);
    % 'SymbolMapping','Custom',...
    %'CustomSymbolMapping',vect,...
qamDemod = comm.RectangularQAMDemodulator(...
    'ModulationOrder',nsymconst, ...
    'BitOutput', true, ...
    'VarianceSource','Property',...
    'Variance', variance,...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'NormalizationMethod','Average power');
datademapping=step(qamDemod, dataaffft);
datademapping1=reshape(datademapping,bps,[]);

hard_decisionnrotncan =datademapping1<0;
hard_decisionnrotncan1=reshape(hard_decisionnrotncan,[],1);

%Temp=reshape(hard_decisionnrotncan1,(2*NFFT+2),[]);
%real_hard_decisionnrotncan1=Temp(2:bps*NFFT+1,:);

To_compare=hard_decisionnrotncan1(:);
c=xor(data,To_compare);


nberror=nnz(c);
BER= [BER nberror/length(data)];
SNRplot=[SNRplot SNR];
end
%k= log2(nsymconst);
%Eb_N0_dB=SNRplot-10*log10(k);
%theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));
figure ;
semilogy(SNRplot, BER, 'b-*','LineWidth',2);
%hold on;
%semilogy(Eb_N0_dB, theoryBer, 'm-o','LineWidth',2);
title ('Bit Error Rate vs SNR BPSK OFDM')
grid on
legend('simulation','theory');
toc
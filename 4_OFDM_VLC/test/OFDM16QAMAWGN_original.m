tic
clear all;
close all;
bps=4;
Weather_impact=1;
nsymconst=2^bps;
NFFT=256;
GI=1/32; %correspond au CP
nbsymbol=1000;
data= randi ([0 1],NFFT*nbsymbol*bps,1);
%vect=[8 9 13 12 10 11 15 14 2 3 7 6 0 1 5 4];
Qammodulator=comm.RectangularQAMModulator('ModulationOrder', nsymconst, ...
    'BitInput',true,...
    'NormalizationMethod', 'Average power'); 
      %  'SymbolMapping','Custom',...
    %'CustomSymbolMapping',vect,...
dataoutmapping=step(Qammodulator,data);
datanbcnbs=reshape(dataoutmapping,NFFT,nbsymbol);
dataOFDM=ifft(datanbcnbs,NFFT);
dataifft=dataOFDM([(end-(NFFT*GI-1)):end 1:end],:);
zerosOFDMreals=zeros(2*size(dataifft,1),nbsymbol);
for i=1:nbsymbol
    datareal=real(dataifft(:,i));
    dataim=imag(dataifft(:,i));
    zerosOFDMreals(1:2:2*size(dataifft,1),i)=datareal;
    zerosOFDMreals(2:2:2*size(dataifft,1),i)=dataim;
    %zerosUFMCreal(:,i)=cat(1,real(UFMCcomplex),imag(UFMCcomplex));  
end
ZerosUFMCrealspositive= zerosOFDMreals;
ZerosUFMCrealsnegative= zerosOFDMreals;
DCB=7;%dB
zeroDCOFDM=zeros(2*size(dataifft,1),nbsymbol);
        for symIdx=1:nbsymbol
            %txSigt= zerosUFMCreals(:,symIdx);
         p1 = sum(abs(zerosOFDMreals(:,symIdx)).^2)/length(zerosOFDMreals(:,symIdx));%calcul de la puissance moyenne
        % noise generation
        %DCBiais=sqrt(10^(DCB/10)-1)*sqrt(p1);
        DCBiais=0.5;
       %DCBiais=0.1;
        zeroDCOFDM(:,symIdx)=(DCBiais+zerosOFDMreals(:,symIdx))*Weather_impact;
        end


dataifft1=zeroDCOFDM(:);
BER=[];
SNRplot=[];
for SNR=0:1:40
    
p = sum(dataifft1.*conj(dataifft1))/length(dataifft1);
snr = 10^(SNR/10);%valeur linéaire du rapport signal à bruit
desv = sqrt(p/snr)/sqrt(2);
n = randn(size(dataifft1));% %+ 1i*randn(size(dataifft1));
n = n*desv;
signalnoisy= dataifft1 + n; 


signalnoisyy=reshape(signalnoisy,2*(NFFT+GI*NFFT),[]);
    zeroDCOFDMrec=zeros(size(signalnoisyy,1),nbsymbol);
        for symIdx=1:nbsymbol
            %txSigt= UFMCcomplexr(:,symIdx);
         p2 = sum(abs(signalnoisyy(:,symIdx)).^2)/length(signalnoisyy(:,symIdx));%calcul de la puissance moyenne
        % noise generation
       %DCBiais2=0.1;
        %DCBiais2=sqrt(10^(DCB/10)-1)*sqrt(p2);
        DCBiais2=0.5;
        zeroDCOFDMrec(:,symIdx)=-DCBiais2+signalnoisyy(:,symIdx);
        end
        %for i=1:numUFMCSymbols
            %UFMCdata=UFMCcomplexr(:,i);
    %UFMCdatawithoutht=UFMCdata(:,1:2*(numFFT+filterLen-1));
    %rxSig=rxSig1(1:length(txSig));
    zerosOFDMrealsR=zeros(size(dataifft,1),nbsymbol);
    datapostOQAMr=zeros(2*(NFFT+GI*NFFT),nbsymbol);
   for i=1:nbsymbol
                      datapostOQAMr(1:2:2*(NFFT+GI*NFFT),i)= real(zeroDCOFDMrec(1:2:2*(NFFT+GI*NFFT),i));
                datapostOQAMr(2:2:2*(NFFT+GI*NFFT),i)= 1i*real(zeroDCOFDMrec(2:2:2*(NFFT+GI*NFFT),i));
                zerosOFDMrealsR(:,i)=datapostOQAMr(1:2:2*(NFFT+GI*NFFT),i)+datapostOQAMr(2:2:2*(NFFT+GI*NFFT),i);
         
   end



zerosOFDMrealsR1=zerosOFDMrealsR(NFFT*GI+1:NFFT + NFFT*GI,:);
datafft=fft(zerosOFDMrealsR1,NFFT);  
dataaffft=datafft(:);
variance=10^(-SNR/10);
    % 'SymbolMapping','Custom',...
    %'CustomSymbolMapping',vect,...
% qamDemod = comm.RectangularQAMDemodulator(...
%     'ModulationOrder',nsymconst, ...
%     'BitOutput', true, ...
%     'VarianceSource','Property',...
%     'Variance', variance,...
%     'DecisionMethod','Approximate log-likelihood ratio',...
%     'NormalizationMethod','Average power');
qamDemod = comm.RectangularQAMDemodulator(...
    'ModulationOrder',nsymconst, ...
    'BitOutput', true, ...
    'DecisionMethod','Approximate log-likelihood ratio',...
    'NormalizationMethod','Average power');
datademapping=step(qamDemod, dataaffft);
datademapping1=reshape(datademapping,bps,[]);
hard_decisionnrotncan =datademapping1<0;
hard_decisionnrotncan1=reshape(hard_decisionnrotncan,[],1);
c=xor(data,hard_decisionnrotncan1);
nberror=nnz(c);
BER= [BER nberror/length(data)];
SNRplot=[SNRplot SNR];
end
%k= log2(nsymconst);
%Eb_N0_dB=SNRplot-10*log10(k);
%theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));
%figure (1);
sn=SNRplot;
be=BER;
semilogy(sn, be, 'b-*','LineWidth',2);
%hold on;
%semilogy(Eb_N0_dB, theoryBer, 'm-o','LineWidth',2);
title ('Bit Error Rate vs SNR 16-QAM')
grid on
legend('simulation','theory');
hold on;
toc
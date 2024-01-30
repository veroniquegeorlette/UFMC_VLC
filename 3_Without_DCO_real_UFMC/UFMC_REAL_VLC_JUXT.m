%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UFMC modulation implementation in a lite version of DVB-T2 system             %                                                            %
%                                                                                %
% Description : UFMC modulation (1K mode), 16-QAM mod, VLC chan,    %                                                     %                                                                                                                                 %                                                               %
% Zero Forcing equalization                                                      %
% Author : Anne-Carole HONFOGA & Véronique Georlette                                                  %                                                                      %
% Date   : Janvier 2021                                                          %                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
numFFT = 256;%1024;        % number of FFT points
subbandSize =12;%12 %66;%14;    % must be > 1
numSubbands =21; %21 %14;%64;    % numSubbands*subbandSize <= numFFT
subbandOffset = 3;%2;%50;%64 % numFFT/2-subbandSize*numSubbands/2 for band center
numUFMCSymbols=10;
% Dolph-Chebyshev window design parameters
filterLen = 18;%36;%32%172;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB

bitsPerSubCarrier = 1;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
%SNR = 10;              % SNR in dB
% Design window with specified attenuation
prototypeFilter = chebwin(filterLen, slobeAtten);

qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%  Initialize arrays
inpData = zeros(numSubbands,bitsPerSubCarrier*subbandSize );
inputData=zeros(bitsPerSubCarrier*subbandSize*numSubbands,numUFMCSymbols);
%txSig = complex(zeros(numFFT+filterLen-1, 1));
UFMCcomplex = complex(zeros(numFFT+filterLen-1, numUFMCSymbols));
%hFig = figure(1);

bitsInSym=randi([0 1], bitsPerSubCarrier*subbandSize*numSubbands, numUFMCSymbols);


%[ht,delay]=rayleigh_channel_profile('TU6',numFFT);
for symIdx=1:numUFMCSymbols
%  Loop over each subband
    txSig=complex(zeros(numFFT+filterLen-1,1));
    dataeachband=complex(zeros(numFFT+filterLen-1,numSubbands));
    
    
    bitsInSym_ = bitsInSym(:,symIdx);
    symbolsIn = qamMapper(bitsInSym_);
    bitsInSym_resh=reshape(symbolsIn,subbandSize,numSubbands); %subband size of lines

    for bandIdx = 1:numSubbands
   
        
        
         bitsIn = bitsInSym_resh(:,bandIdx);
        bitsIn_=bitsIn';
        inpData(bandIdx,:) = bitsIn_; % log bits for comparison

        % Pack subband data into an OFDM symbol
        offset = subbandOffset+(bandIdx-1)*subbandSize;
        symbolsInOFDM = [zeros(offset,1); bitsIn; ...
                         zeros(numFFT-offset-subbandSize, 1)];
 %this is a real signal
        ifftOut = ifft(ifftshift(symbolsInOFDM));

        % Filter for each subband is shifted in frequency
        bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                     ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );
% This is not a real signal
        filterOut = conv(bandFilter,ifftOut);

        % Plot power spectral density (PSD) per subband
        [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                              numFFT*2, 1, 'centered');
        plot(f,10*log10(psd));

        % Sum the filtered subband responses to form the aggregate transmit
        % signal
        txSig = txSig + filterOut;
        dataeachband(:,bandIdx)= filterOut;
    end
    UFMCcomplex(:,symIdx)=txSig;
    inputData(:,symIdx)=bitsInSym_(:);
end
%% real UFMC signal generation
%zerosUFMCreal=zeros(2*size(UFMCcomplex,1),numUFMCSymbols);
zerosUFMCreals=zeros(2*size(UFMCcomplex,1),numUFMCSymbols);
for i=1:numUFMCSymbols
    datareal=real(UFMCcomplex(:,i));
    dataim=imag(UFMCcomplex(:,i));
    zerosUFMCreals(1:2:2*size(UFMCcomplex,1),i)=datareal;
    zerosUFMCreals(2:2:2*size(UFMCcomplex,1),i)=dataim;
    %zerosUFMCreal(:,i)=cat(1,real(UFMCcomplex),imag(UFMCcomplex));  
end
ZerosUFMCrealspositive= zerosUFMCreals;
ZerosUFMCrealsnegative= zerosUFMCreals;
%DCB=13;%dB
zeroDCUFMC=zeros(2*size(UFMCcomplex,1),numUFMCSymbols);
        for symIdx=1:numUFMCSymbols
            %txSigt= zerosUFMCreals(:,symIdx);
         p = sum(abs(zerosUFMCreals(:,symIdx)).^2)/length(zerosUFMCreals(:,symIdx));%calcul de la puissance moyenne
        % noise generation
        DCBiais=0.5;
        %DCBiais=sqrt(10^(DCB/10)-1)*sqrt(p);
        %DCBiais=1;
        zeroDCUFMC(:,symIdx)=zerosUFMCreals(:,symIdx);%+DCBiais;
        end
% for i=1:numUFMCSymbols
%     for j=1:size(ZerosUFMCrealspositive,1)
%         if ZerosUFMCrealspositive(j,i)<0
%             ZerosUFMCrealspositive(j,i)=0;
%         else
%             ZerosUFMCrealspositive(j,i)=ZerosUFMCrealspositive(j,i);
%         end
%     end
% end 
% 
% for i=1:numUFMCSymbols
%     for j=1:size(ZerosUFMCrealsnegative,1)
%         if ZerosUFMCrealsnegative(j,i)<0
%             ZerosUFMCrealsnegative(j,i)=0;
%         else
%             ZerosUFMCrealsnegative(j,i)=ZerosUFMCrealsnegative(j,i);
%         end
%     end
% end 

% for i=1:numUFMCSymbols
%     for j=1:size(zerosUFMCreals,1)
% Positivesvalues=zerosUFMCreals(i,j)
% 
% S01=sousregion1norm(sousregions1binairebi(:,i)==0);
%                 S11=sousregion1norm(sousregions1binairebi(:,i)==1);
                
%set(hFig, 'Position', figposition([20 50 25 30]));
%hold off;

% Plot power spectral density (PSD) per subband
%figure (2)
%[psd,f] = periodogram(dataeachband(:,10), rectwin(length(filterOut)), ...
                             % numFFT*2, 1, 'centered');
%plot(f,10*log10(psd),'b');
%grid on
%axis([-0.5 0.5 -100 20]);
            %% Plot power spectral density (PSD) for symbol
[psd,f] = periodogram(zerosUFMCreals(1:(numFFT+filterLen-1)/2,1), rectwin(size(zerosUFMCreals(1:(numFFT+filterLen-1)/2,1),1)), ...
                              numFFT*4, 1,'centered');
Nbporteuse=subbandSize*numSubbands;
figure (3)
plot(f,10*log10(psd),'b');
grid on
%axis([-0.5 0.5 -100 20]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['UFMC, ' num2str(Nbporteuse) ' Subcarriers'])

% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txSig);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);
%BERDEM=[];
%BERREAL=[];
%NBrealisation=100;
%for nb=1:NBrealisation
%ht=rayleigh_channel_profile('OneTap',numFFT);
%hF = fftshift(fft(real(ht),numFFT,2));
BER1=[];
SNRplot=[];
%[ht,delay]=rayleigh_channel_profile('TU6',numFFT);
    for SNR=0:1:30
    %txSig1=conv( zerosUFMCreals,ht);
    %rxSig = awgn(txSig, SNR, 'measured');
   %Signalt1=zeros(numUFMCSymbols,2*(numFFT+filterLen-1)+length(ht)-1);
   UFMCcomplextot= zeroDCUFMC.';
     %UFMCcomplextot=zerosUFMCreals.';
    %for symIdx=1:numUFMCSymbols
          %Signalt1(symIdx,:)=conv(UFMCcomplextot(symIdx,:),real(ht));
    %end
        %UFMCcomplexafch=Signalt1.';
        UFMCcomplexafch=UFMCcomplextot.';
        nbcarrier = subbandSize*numSubbands;
        TU=112e-06;% durée du symbole
         %Signal bandwidth calcul de la bande passante du signal
        sBW = nbcarrier/TU;% elle représente la fréquence d'échantillonnage% nbporteuse doit représenter le nombre de sous porteuse de donnée
        % Noise bandwidth
        nBW = numFFT/TU;% numFFT nombre de sous porteuse totale, le calcul de la bande passante du bruit inclu le nombre de sous porteuse totale
        UFMCcomplexr = zeros(2*(numFFT+filterLen-1), numUFMCSymbols);
        for symIdx=1:numUFMCSymbols
            txSigt= UFMCcomplexafch(:,symIdx);
            %p = sum(abs(zerosUFMCreals(:,symIdx)).^2)/length(zerosUFMCreals(:,symIdx));
         p = sum(abs(UFMCcomplexafch(:,symIdx)).^2)/length(UFMCcomplexafch(:,symIdx));%calcul de la puissance moyenne
        % noise generation
         snr = 10^(SNR/10);%valeur linéaire du rapport signal à bruit
         No=p/snr;
           desv = sqrt(No)/sqrt(sBW/nBW);
        n = desv*(randn(size(txSigt)));% + 1i*randn(size(txSigt)));
        %n = n*desv;
        rxSig=txSigt+n;
        UFMCcomplexr(:,symIdx)=rxSig;
        end
        
        zeroDCUFMCrec=zeros(2*size(UFMCcomplex,1),numUFMCSymbols);
        for symIdx=1:numUFMCSymbols
            %txSigt= UFMCcomplexr(:,symIdx);
         p = sum(abs( UFMCcomplexr(:,symIdx)).^2)/length( UFMCcomplexr(:,symIdx));%calcul de la puissance moyenne
        % noise generation
        %DCBiais2=1;
        %DCBiais2=sqrt(10^(DCB/10)-1)*sqrt(p);
        %DCBiais2=2;
        zeroDCUFMCrec(:,symIdx)=UFMCcomplexr(:,symIdx);%-DCBiais2;
        end
        %for i=1:numUFMCSymbols
            %UFMCdata=UFMCcomplexr(:,i);
    %UFMCdatawithoutht=UFMCdata(:,1:2*(numFFT+filterLen-1));
    %rxSig=rxSig1(1:length(txSig));
    zerosUFMCrealsR=zeros(size(UFMCcomplex,1),numUFMCSymbols);
    datapostOQAMr=zeros(2*(numFFT+filterLen-1),numUFMCSymbols);
   for i=1:numUFMCSymbols
                      datapostOQAMr(1:2:2*(numFFT+filterLen-1),i)= real(zeroDCUFMCrec(1:2:2*(numFFT+filterLen-1),i));
                datapostOQAMr(2:2:2*(numFFT+filterLen-1),i)= 1i*real(zeroDCUFMCrec(2:2:2*(numFFT+filterLen-1),i));
                zerosUFMCrealsR(:,i)=datapostOQAMr(1:2:2*(numFFT+filterLen-1),i)+datapostOQAMr(2:2:2*(numFFT+filterLen-1),i);
         
   end
    
% for i=1:numUFMCSymbols
%     dataR=UFMCcomplexr(:,i);
%     datad=zerosUFMCrealsR(:,i)
%     dataht=dataR(:,2*numFFT+1:size(UFMCcomplexr,1));
%     zerosUFMCcomplexreceived(:,i)=cat(1,datad,dataht);  
% end
    rxBits=zeros(bitsPerSubCarrier*subbandSize*numSubbands,numUFMCSymbols);
     for symIdx=1:numUFMCSymbols
         yRx=zerosUFMCrealsR(:,(symIdx));

        % Pad receive vector to twice the FFT Length (note use of txSig as input)
        %   No windowing or additional filtering adopted
        yRxPadded = [yRx; zeros(2*numFFT-numel(txSig),1)];

        % Perform FFT and downsample by 2
        RxSymbols2x = fftshift(fft(yRxPadded));
        RxSymbols = RxSymbols2x(1:2:end);
        %channel equalization
       % RxSymbolsafeq=RxSymbols.'./real(hF);
        % Select data subcarriers
        RxSymbolsafeq=RxSymbols.';
        %dataRxSymbols = RxSymbolsafeq(subbandOffset+(1:numSubbands*subbandSize));
        dataRxSymbols = RxSymbolsafeq(subbandOffset+(1:numSubbands*subbandSize));
        % Plot received symbols constellation
        %%constDiagRx = comm.ConstellationDiagram('ShowReferenceConstellation', ...
            %%false, 'Position', figposition([20 15 25 30]), ...
           %% 'Title', 'UFMC Pre-Equalization Symbols', ...
           %% 'Name', 'UFMC Reception', ...
            %%'XLimits', [-150 150], 'YLimits', [-150 150]);
        %%constDiagRx(dataRxSymbols);

        % Use zero-forcing equalizer after OFDM demodulation
        rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
               zeros(numFFT-filterLen,1)];
        prototypeFilterFreq = fftshift(fft(rxf));
        prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));
        % Equalize the channel
        
        % Equalize per subband - undo the filter distortion
        dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
        EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
        EqualizedRxSymbols = EqualizedRxSymbolsMat(:);

        % Plot equalized symbols constellation
    %     constDiagEq = comm.ConstellationDiagram('ShowReferenceConstellation', ...
    %         false, 'Position', figposition([46 15 25 30]), ...
    %         'Title', 'UFMC Equalized Symbols', ...
    %         'Name', 'UFMC Equalization');
    %     constDiagEq(EqualizedRxSymbols);
    %             'SymbolMapping','Custom',...
    %      'CustomSymbolMapping',vect,...
        % Demapping and BER computation
        
        qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
            2^bitsPerSubCarrier, 'BitOutput', true, ...
            'NormalizationMethod', 'Average power', 'Variance',10^(-SNR/10));
        BER = comm.ErrorRate;
    rxBits(:,(symIdx))=qamDemod(EqualizedRxSymbols);
     end
    % Perform hard decision and measure errors
    %rxBits = qamDemod(EqualizedRxSymbols);
    ber = BER(inputData(:), rxBits(:));
    BER1=[BER1 ber(1)];
    SNRplot=[SNRplot,SNR];
    %disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
      %  num2str(SNR) ' dB']);
    end
 %BERDEM=[BERDEM;BER1];
 %BERREAL=mean(BERDEM);
% end
figure (4);
semilogy(SNRplot, BER1,'b-*','LineWidth',2);
grid on;
title ('Bit Error Rate vs SNR BPSK UFMC')
hold on;
%figure (5)
%semilogy(SNRplot, BER1);
%title ('Bit Error Rate vs SNR 16-QAM')
%semilogy(SNRplot, BER1)
%Restore RNG state
%rng(s);
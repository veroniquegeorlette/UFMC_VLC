%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UFMC modulation implementation in a lite version of DVB-T2 system             %                                                            %
%                                                                                %
% Description : UFMC modulation (1K mode), 16-QAM mod, VLC chan,    %                                                     %                                                                                                                                 %                                                               %
% Zero Forcing equalization                                                      %
% Author : Anne-Carole HONFOGA                                                   %                                                                      %
% Date   : Janvier 2021                                                          %                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
numFFT = 256;%1024;        % number of FFT points
subbandSize =12;%12 %66;%14;    % must be > 1
numSubbands =21; %21 %14;%64;    % numSubbands*subbandSize <= numFFT
subbandOffset = 3;%2;%50;%64 % numFFT/2-subbandSize*numSubbands/2 for band center
numUFMCSymbols=1;
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
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
inputData=zeros(bitsPerSubCarrier*subbandSize*numSubbands,numUFMCSymbols);
%txSig = complex(zeros(numFFT+filterLen-1, 1));
UFMCcomplex = complex(zeros(2+2*numFFT+filterLen-1, numUFMCSymbols));

bitsInSym=randi([0 1], bitsPerSubCarrier*subbandSize*numSubbands, numUFMCSymbols);


%[ht,delay]=rayleigh_channel_profile('TU6',numFFT);
for symIdx=1:numUFMCSymbols
%  Loop over each subband
    txSig=complex(zeros(2+2*numFFT+filterLen-1,1));
    dataeachband=complex(zeros(2+2*numFFT+filterLen-1,numSubbands));
    bitsInSym_=bitsInSym(:,symIdx);
    bitsInSym_resh=reshape(bitsInSym_,subbandSize,[]);
    for bandIdx = 1:numSubbands
   bitsIn = bitsInSym_resh(:,numUFMCSymbols);
   symbolsIn = qamMapper(bitsIn);
   inpData(:,bandIdx) = bitsIn; % log bits for comparison

        % Pack subband data into an OFDM symbol
        offset = subbandOffset+(bandIdx-1)*subbandSize;
        symbolsInUFMC = [zeros(offset,1); symbolsIn; ...
                         zeros(numFFT-offset-subbandSize, 1)];
                                        
% Need to put Hermitian Symmetry here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
symbolsInUFMC_to_HS=[0,symbolsInUFMC.',0,fliplr(conj(symbolsInUFMC.'))];
temp=symbolsInUFMC_to_HS.';
        ifftOut = ifft(ifftshift(temp));

        % Filter for each subband is shifted in frequency
        bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                     ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );
        filterOut = conv(bandFilter,ifftOut);
filterOut_RE=real(filterOut);
        % Plot power spectral density (PSD) per subband
        [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                              numFFT*2, 1, 'centered');
        plot(f,10*log10(psd));

        % Sum the filtered subband responses to form the aggregate transmit
        % signal
        txSig = txSig + filterOut_RE;
        dataeachband(:,bandIdx)= filterOut_RE;
        data_tot(:,bandIdx)= symbolsIn;
    end

    UFMCcomplex(:,symIdx)=txSig;
    inputData(:,symIdx)=inpData(:);
end



%% Plot power spectral density (PSD) for symbol
%Need to check if this line is ok !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[psd,f] = periodogram(UFMCcomplex(1:(numFFT+filterLen-1)/2,1), rectwin(size(UFMCcomplex(1:(numFFT+filterLen-1)/2,1),1)), ...
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
[~,~,paprUFMC] = PAPR(UFMCcomplex);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);


BER1=[];
SNRplot=[];
%[ht,delay]=rayleigh_channel_profile('TU6',numFFT);
for SNR=0:1:30        
   %UFMCcomplextot= UFMCcomplex.';
   UFMCcomplexafch=UFMCcomplex;
        nbcarrier = subbandSize*numSubbands;
        TU=112e-06;% durée du symbole
         %Signal bandwidth calcul de la bande passante du signal
        sBW = nbcarrier/TU;% elle représente la fréquence d'échantillonnage% nbporteuse doit représenter le nombre de sous porteuse de donnée
        % Noise bandwidth
        nBW = numFFT/TU;% numFFT nombre de sous porteuse totale, le calcul de la bande passante du bruit inclu le nombre de sous porteuse totale
        %UFMCcomplexr = zeros(2*(numFFT+filterLen-1), numUFMCSymbols);
        UFMCcomplexr = zeros(2+2*numFFT+filterLen-1, numUFMCSymbols);

        
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
        %on rajoute le bruit
        rxSig=txSigt;%+n;
        UFMCcomplexr(:,symIdx)=rxSig;
        end
        
              
   % empty element
    rxBits=zeros(bitsPerSubCarrier*subbandSize*numSubbands,numUFMCSymbols);
    
    
     for symIdx=1:numUFMCSymbols
         yRx=UFMCcomplexr(:,(symIdx));

        % Pad receive vector to twice the FFT Length (note use of txSig as input)
        %   No windowing or additional filtering adopted
        yRxPadded = [yRx; zeros(4*numFFT-numel(txSig),1)];

        % Perform FFT and downsample by 2
        RxSymbols2x = fftshift(fft(yRxPadded));
        Temp=RxSymbols2x(2:2*(numFFT+1),:);
        RxSymbols = Temp(1:2:end);
        
        %Temp=RxSymbols(2:numFFT+1,:)';
        
        
      
        % Select data subcarriers
        
        %dataRxSymbols = RxSymbolsafeq(subbandOffset+(1:numSubbands*subbandSize));
        dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));
        

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
        %Temp=EqualizedRxSymbols(2:numFFT+1,:);
        
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

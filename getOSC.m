%OTC example
%edit 25/9/2014
%Copyright (c) Dino Dvorak neuraldino@gmail.com
%
%code converts time series data into time-frequency representation,
%normalizes it and finds peaks above predefined power threshold. frequency
%specific peaks are then used as trigger points for time domain averaging
%of the raw signal. resulting signal is the OTCG - OTC comodulogram - which
%shows modulatory frequencies, modulated frequencies and phase of their
%coupling
%Citation:
%J Neurosci Methods. 2014 Mar 30;225:42-56. doi: 10.1016/j.jneumeth.2014.01.002. Epub 2014 Jan 19.
%Toward a proper estimation of phase-amplitude coupling in neural oscillations.
%Dvorak D, Fenton AA
clear all; close all; clc;

%load file
load('lfp_rat_hpc.mat');
%normalize
eegData = eegData / std(eegData);

%bands for finding oscillations
bandsOSC = 10:1:130;
NbandsOSC = length(bandsOSC);

%bands for finding modulation
bandsOTC = 20:120;
NbandsOTC = length(bandsOTC);

powerTh = 2; %power threshold (S.D.)

winLen = 0.6*eegFS; %OTC window

winProc = 10*eegFS; %processing window

tic

%create set of Morlet wavelets
wFactor = 8; %scale factor
wavelets = cell(1,NbandsOSC);
for bI = 1:length(bandsOSC)
    f = bandsOSC(bI);
    sigmaF = f/wFactor;    %practical setting for spectral bandwidth (Tallon-Baudry)
    sigmaT = 1/(sigmaF*pi*2); %wavelet duration
    t = -4*sigmaT*eegFS:4*sigmaT*eegFS;
    t = t/eegFS;
    if rem(length(t),2) == 0; t = t(1:end-1); end
    S1 = exp((-1*(t.^2)) / (2*sigmaT^2));
    S2 = exp(2*1i*pi*f*t);
    A = (sigmaT * sqrt(pi))^(-0.5); %normalization
    psi = A*S1.*S2;
    wavelets{bI} = psi;
end

%number of processing windows
Nwin = floor(length(eegData)/winProc);

%store samples and frequencies of high power oscillations
osc = [];
for wI = 1:Nwin

    stW = (wI-1)*winProc + 1;
    edW = stW + winProc - 1;

    eeg = eegData(stW:edW);

    pwrN = zeros(NbandsOSC,winProc); %z-score power

    %compute nomalized power across bands
    for bandI = 1:NbandsOSC
        psi = wavelets{bandI};
        %convolution
        c = conv(eeg,psi);
        N = round((length(psi)-1)/2);
        c = c(N:length(c)-N); 
        if length(c) > size(eeg,2); c = c(1:size(eeg,2)); end;

        %power
        power = (abs(c)).^2;
        %normalized power
        pwrN(bandI,1:length(power)) = power / (std(c))^2;
    end

    %find local maxima in normalized power profile
    bw = imregionalmax(pwrN);
    [y,x] = find(bw==1);
    pN = pwrN(bw);

    %extract powers higher than threshold not on the border
    k = pN > powerTh & x > 1 & x < size(pwrN,2) & y > 1 & y < size(pwrN,1);
    x = x(k);
    y = y(k);

    %convert to samples and frequencies
    samples = x + stW - 1;
    freq = bandsOSC(y);

    osc = cat(1,osc,cat(2,samples,freq'));
end

%for each band, extract relevant samples and average raw LFP around them
otc = zeros(NbandsOTC,winLen*2);
for bI = 1:NbandsOTC
    f = bandsOTC(bI);
    stB = bandsOTC(bI)-f/wFactor; %bandwidth proportional to freuquency
    edB = bandsOTC(bI)+f/wFactor;
    
    %find samples at specific frequency
    k = osc(:,2) >= stB & osc(:,2) < edB;
    samples = osc(k,1);
    samples = samples(samples > winLen & samples < length(eegData)-winLen);
    
    aver = zeros(1,winLen*2);
    
    for sI = 1:length(samples)
        eeg = eegData(samples(sI)-winLen+1:samples(sI)+winLen);
        aver = aver + eeg;
    end
    
    otc(bI,1:length(aver)) = aver;
end

%compute modulation strength profile
modS = zeros(NbandsOTC,1);
for bI = 1:NbandsOTC
    o = otc(bI,:);
    m = max(o) - min(o);
    modS(bI) = m;
end
modS = smooth(modS,5); %smooth

toc
    
figure
subplot(1,5,1:4); hold on;
imagesc(otc)
plot([winLen winLen],[1,NbandsOTC],'w','LineWidth',2);
otc45 = otc(bandsOTC == 45,:); otc45 = otc45/max(otc45); otc45 = otc45*5;
otc70 = otc(bandsOTC == 70,:); otc70 = otc70/max(otc70); otc70 = otc70*5;
plot(otc45+(45-bandsOTC(1)),'k');
plot(otc70+(70-bandsOTC(1)),'k');
axis xy tight
set(gca,'XTick',[1,winLen,winLen*2])
set(gca,'XTickLabel',{-1*winLen/eegFS,0,winLen/eegFS})
set(gca,'YTick',1:10:length(bandsOTC))
set(gca,'YTickLabel',{bandsOTC(1:10:end)})
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
title('OTC comodulogram');

subplot(1,5,5); hold on;
plot(modS,1:length(modS),'k');
axis xy tight
set(gca,'XTickLabel',[])
set(gca,'YTick',1:10:length(bandsOTC))
set(gca,'YTickLabel',{bandsOTC(1:10:end)})
xlabel('Mod. strength');


 
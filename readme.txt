Copyright (C) Dino Dvorak 2014 neuraldino@gmail.com
content:
getOSC.m - code for computing OTC comodulogram as per Dvorak & Fenton JNMethods 2014 paper
lfp_rat_hpc.mat - sample file with LFP recording from rat CA1 part of hippocampus
otc.pes - results of the calculation - OTC comodulogram

The code computes normalized time-frequency power profile using Morlet wavelets. It then searches for local maxima (peaks) and filters only these with large enough power (e.g. > 2 S.D.). These maxima are then selected in band-specific manner and used as trigger points for time-domain averaging of windowed LFP signal. If phase-amplitude coupling is present, peaks of detected frequency-specific oscillations (gamma, 20-100Hz) would be systematically locked to a slow frequency oscillation (theta, 8Hz) and therefore destructive averaging would remove all non-phase locked signal while the constructive averaging would amplify the signal, phase locked to peaks of detected oscillations. The resulting signal is the modulatory signal (theta). Its phase in the middle of the window (time 0) corresponds to the modulatory phase (location of gamma peak). Amplitude of the modulatory signal (min-max) corresponds to the relative strength of coupling and is shown on the right side of the attached figure. Clear slow- and fast-gamma peaks around 45 and 70 Hz respectively are clearly visible. 

Citation of work:
J Neurosci Methods. 2014 Mar 30;225:42-56. doi: 10.1016/j.jneumeth.2014.01.002. Epub 2014 Jan 19.
Toward a proper estimation of phase-amplitude coupling in neural oscillations.
Dvorak D, Fenton AA
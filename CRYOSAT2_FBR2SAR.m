% CRYOSAT2_FBR2SAR.m : is a MATLAB program that perfroms SAR processing of
% CryoSat-2 SAR altimeter data. Please see Wingham et al., 2006
% (included in this software bundle), for a detailed description of the 
% processing done to the received signal to form the L1B echo. 
% This code ingests the bursts from FBR data, forms the beams from the 
% pulses within each burst, corrects the beams for the altitude rate, and 
% applies a slant-range correction to align the beams in range. 
% If these beams are summed across a burst, over ocean the sum will 
% approximate an L1B ocean echo. However, it won't be exact, as this code
% does not perform stacking of the beams.
%
% Copyright (C) 2013 Natalia Galin.
% Questions/Comments/Suggestion/Corretions? 
% Please do get in touch: natalia.galin@gmail.com
% (A doc describing each useful line of this code is in the pineline, and 
% will be released shortly. So stay tuned if you're intrested in the 
% nitty gritty.)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%-----------------------------------------------------------------------%%
%% INPUTS REQUIRED:
%  COMPLEX_DATA_I - [RECNO x BURSTNO x PULSEPERBURST x SAMPLENO] matrix where:
%                   RECNO: is the number of 'records' (the FBR data is stored in batches of
%                          20 bursts, each set of 20 bursts is called a 'record')
%                   BURSTNO: is '20'. 
%                   PULSEPERBURST: is '64', the number of pulses per burst (see Wingham et al., 2006)
%                   SAMPLENO: is 128, the number of samples per burst
%
%  COMPLEX_DATA_Q - (as 'COMPLEX_DATA_I')
%
%  SAT_VEL        - [RECNO], vector of satellite velocities, updated once every 20
%                   bursts, units must be [m/s]
%
%  ALT_RATE       - [RECNO], vector of altitude rate of the satellite,
%                   updated once every 20 bursts, units must be [m/s]
%
%  ALTITUDE       - [RECNO], vector of the satellite altitude, updated once
%                   every 20 bursts, units must be [m]
% 
%% DEPENDENCIES:
%  FUNC_OCOG_SARv1.m - function that calculates an empirical estimate for
%                      the 'point of closest approach' based on the OCOG 
%                      method, see Wingham et al., 1986. 
%
%  interpResampFFTSAR.m - function that resamples the waveform to center 
%                         its 'point of closest approach' in the center
%                         of the range window.
%
%% OPTIONS/FLAGS:
%  RETRACK_FLAG - set to '1' if you'd like to have each pulse 'retracked',
%                 i.e resampled such that each pulse's empirically determined 
%                 'point of first return' is aligned to the centre point of 
%                 the range window. As the range window is increased by the
%                 code from 128 to 256 bins. The increase necessary in order
%                 to preserve the frequency content of the pulse power which 
%                 is formed by taking the squared modulus of 'I+sqrt(-1)*Q'. 
%               - set to '0' if you don't want to have each pulse 'retracked'.   
%
%  HAMWIN_FLAG  - set to '1' if you would like to have a Hamming window
%                 applied across the 64 pulses within each burst. 
%               - set to '0' if you don't want to have the Hamming window
%                 applied.
%                 Hamming Window Definition:
%                 n = [0:64]; 
%                 hamWin = 0.54 - 0.46*cos(2*pi*n/(N-1));
%
%% OUTPUTS: 
%  BURSTPOWER - [TOTALBURSTS x BEAMS x UPSAMPLENO] matrix where:
%               TOTALBURSTS: is equal to the total number of bursts, i.e. 
%                            TOTALBURSTS = RECNO x BURSTNO
%               BEAMS: is 64, each pulse is now changed to a beam
%               UPSAMPLENO: is 128*2 = 256, to preserve the frequency
%               content of the beam power.
%%-----------------------------------------------------------------------%%                

function [BURSTPOWER] = CRYOSAT2_FBR2SAR(COMPLEX_DATA_I,COMPLEX_DATA_Q,SAT_VEL,ALT_RATE,ALTITUDE,RETRACK_FLAG,HAMWIN_FLAG)

[RECNO,BURSTNO,PULSEPERBURST,SAMPLENO] = size(COMPLEX_DATA_I);

TOT_RECORD = RECNO;
STR_RECORD = 1;

%DECLARE CONSTANTS:
BLOCK = 20;
ECHO_IN_BURST = 64;
SAMPLE = 128;
TRK_WIN = 60; %[m]
N = 64;
%DEFINE: hamming window
n = [1:N]-1;
hamWin = 0.54 - 0.46*cos(2*pi*n/(N-1));
PRF = 18181.82;
BRI = 11.7e-3; %for SAR mode
c = 299792458; 
lam = c/13.575e9; 
k0 = 2*pi/lam;
R = 6378.1370e3; %radius of Earth in meters

%now some rudimentary error checking to make sure your structures are at least
%in the correct format - can't be sure about the data ;)
if (BURSTNO ~= BLOCK)
    exit;
end
if (PULSEPERBURST ~= ECHO_IN_BURST)
    exit;
end
if (SAMPLENO ~= SAMPLE)
    exit;
end

%SOME INITIALISATION:
idxBurst = 0;
interp_chan1 = zeros(1,ECHO_IN_BURST,SAMPLE);
recordAvePow = zeros(ECHO_IN_BURST,SAMPLE*2);
recordAvePowNorm = zeros(ECHO_IN_BURST,SAMPLE*2);

for recN=STR_RECORD:TOT_RECORD
    SAMPLE = 128;
    recN
    interp_chan1 = zeros(BLOCK,ECHO_IN_BURST,SAMPLE);
    interp_chan1(:,:,:) = squeeze(COMPLEX_DATA_I(recN,:,:,:)) + sqrt(-1)*squeeze(COMPLEX_DATA_Q(recN,:,:,:));
    fft_interp_chan1 = interp_chan1;

    %-------------------------------------------------------------------------%
    %RANGE FFT:
    Ns = SAMPLE;
    clear case4_out1 case4_out2
    delR = TRK_WIN/SAMPLE;
    
    for j=1:BLOCK
        for m=1:ECHO_IN_BURST
            
            fft_interp_chan1_ROT = squeeze(fft_interp_chan1(j,m,:));
            wfm = fftshift(fft_interp_chan1_ROT,1);

            out_wfm = [wfm(1:Ns/2).' zeros(1,Ns) wfm(Ns/2+1:end).'];
            out_wfm(Ns/2+1) = wfm(Ns/2+1)/2;
            out_wfm(Ns+Ns/2+1) = wfm(Ns/2+1)/2;
            out_wfm = fft(out_wfm,Ns*2)/sqrt(Ns*2);
            out_wfm = fftshift(out_wfm);
            
            %case4_out1(j,m,:) = interpResampFFTSAR(out_wfm,1,1,-retrackVals(recN,j)/delR,1/2);
            if (RETRACK_FLAG) 
                retrackChan1_filt(recN,j) = FUNC_OCOG_SARv1(abs(out_wfm))*delR;
                case4_out1(j,m,:) = interpResampFFTSAR(out_wfm,1,1,-retrackChan1_filt(recN,j)/delR,1/2);
%            if (RETRACK_FLAG) 
%                case4_out1(j,m,:) = interpResampFFTSAR(out_wfm,1,1,-retrackChan1_filt(recN,j)/delR,1/2);
            else
                case4_out1(j,m,:) = out_wfm;
            end
                        
        end
    end
    SAMPLE = SAMPLE*2;
    interp_chan1 = case4_out1;
    clear fft_interp_chan1 
    
    %AZIMUTH FFT
    clear shiftTest1

    %DEFINE: 'centering the phase' term:
    t = [1:N]-1;
    tmp = exp(+sqrt(-1)*pi*(N-1)/N*t);
    cPhs = tmp(ones(1,SAMPLE),:).';
    
    %DEFINE: rock angle, which now takes care of the altitude rate
    v = mean(SAT_VEL((recN-1)*BLOCK+1:(recN)*BLOCK));
    PRI = 1/PRF;
    x = v*PRI;
    
    beamInc = pi/(x*k0)/ECHO_IN_BURST*(180/pi);
    h = mean(ALTITUDE((recN-1)*BLOCK+1:(recN)*BLOCK));
    eta = 1 + h/R;
    %calculating the th0, as a function of satellite altitude rate:
    D = (ECHO_IN_BURST/PRF)*v;
    Dt = PRI;
    meanAltRate = mean(ALT_RATE((recN-1)*BLOCK+1:(recN)*BLOCK));
    
    %calculating the corrected velocity, minus altRate:
    th0 = asin(meanAltRate/v);
    thApplied = th0;
    thBeam = th0/(beamInc*pi/180);
    thR = -th0;
    for j=1:BLOCK
        clear shiftTest1 
        %perform a shift of samples:
        for m=1:SAMPLE
            
            if (HAMWIN_FLAG)
                test1 = squeeze(interp_chan1(j,:,m)).*hamWin;
            else
                test1 = squeeze(interp_chan1(j,:,m));
            end
            
            shiftTest1(m,:) = test1.*exp(-2*sqrt(-1)*[0:N-1]*k0*Dt*v*thR);
        end
        
        shiftTest1 = shiftTest1.';
        
        %part 2 of the centering phase term: thR != 0
        cPhs_p2 = exp((N-1)*sqrt(-1)*k0*Dt*v*thR).*ones(N,SAMPLE);
                
        fft_interp_chan1(j,:,:) = cPhs.*cPhs_p2.*fftshift(fft(shiftTest1,N,1),1)/sqrt(N);
    end 
        
    %output:
    burstBeam_chan1 = fft_interp_chan1;
    
    %SLANT RANGE CORRECTING THE ECHOES:
    beamAng = ([1:64]-33).*beamInc;
    
    %slant range correction factor for beams:
    tau = (h*eta/c)*(beamAng*pi/180).^2;
    src = tau*c;
    %initialisation:
    src_fft_interp_chan1 = zeros((BLOCK),ECHO_IN_BURST,SAMPLE);
    src_fft_Pow = zeros(BLOCK,ECHO_IN_BURST,SAMPLE);
    norm_src_fft_Pow = zeros(BLOCK,ECHO_IN_BURST,SAMPLE);
    for i=1:ECHO_IN_BURST
        for j=1:BLOCK
            idxBurst = idxBurst + 1;
            
            tmp = squeeze(burstBeam_chan1(j,i,:));
            tmp = interpResampFFTSAR(tmp.',1,1,src(i)/delR/2,1/2);
            src_fft_interp_chan1(j,i,:) = tmp;
            
            src_fft_Pow(j,i,:) = (abs(src_fft_interp_chan1(j,i,:)).^2);
            
            norm_src_fft_Pow(j,i,:) = (abs(src_fft_interp_chan1(j,i,:)).^2);
        end
        idxBurst = idxBurst - BLOCK;
        
    end
    %OUTPUT:
    BURSTPOWER(idxBurst+1:idxBurst+BLOCK,:,:) = norm_src_fft_Pow;
     
    idxBurst = idxBurst + BLOCK;
end

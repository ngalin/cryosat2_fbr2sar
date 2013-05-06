%CURRENTLY IMPLEMENTED TO GIVE 'EXACTLY' THE SAME VALUES AS THE OCOG
%RETRACKER IN IDL WOULD. HOWEVER, THIS IS NOT ENTIRELY CORRECT, AS IDL
%INDEXING STARTS FROM 0 AND MATLAB FROM 1. NG - 25th Jun, 2011


%CS_OCOG_RETRACKER ported from the IDL version: cs_ocog_retracker.pro
%received from Rob and Marco.

%spec copied from the IDL version:
%--------------------------------------------------------------------------
%--------This function implements the OCOG retracker-----------------------
%mandatory input:
% 1) waveforms:     2D array [n_samples, n_waveforms] of power waveforms
%optional inputs:
% 1) threshold:     percentage of the threshold, see OCOG retracker
% algorithm. set to 0.8 by default
% 2) sampleSize:    size of 1 range bin in meters
%output
% 1) retracker offset: is the offset between the middle of the range window
% and the evaluated surface elevation [metres]
%optional output:
% 1) failure index: indexes of the retracker failures [not implemented in
% MATLAB version]
%
%call example:
%retrack_offset =
%cs_ocog_retracker(waveforms_all,threshold=0.8,sample_size=0.470,failure_index=failure_index_TMP)
% 
%History:
%Created by Marco Fornari, June 2008

function [RTpoint] = FUNC_OCOG_SARv1(wfm)

mid_fraction = 0.8;
sampleSize = 0.47;

[numWfm lenWfm] = size(wfm);

%tunable parameters:
win_centre = lenWfm/2; % LRM/SAR = 64; SARIN = 256
firstbin = 14+1; %all modes
lastbin = lenWfm - (13-1); %LRM/SAR = 115; SARIN = 499

%initialisation:
RTpoint = zeros(1,numWfm);

for i=1:numWfm
    
    %OCOG calculation
    bin_sq = 0;
    wfm_sq_sum = 0;
    wfm_Nsq_sum = 0;
    wfm_qad_sum = 0;
    
    bin = wfm(i,:);
    for j=firstbin:lastbin
        bin_sq = bin(j)*bin(j);
        wfm_sq_sum = wfm_sq_sum + bin_sq;
        wfm_Nsq_sum = wfm_Nsq_sum + (j-1)*bin_sq;
        wfm_qad_sum = wfm_qad_sum + (bin_sq*bin_sq);
    end
    
    %all in units of fractional bins
    cog_posn = wfm_Nsq_sum / wfm_sq_sum;
    OCOG_amp = sqrt(wfm_qad_sum / wfm_sq_sum);
    OCOG_width = (wfm_sq_sum * wfm_sq_sum) / wfm_qad_sum;
    
    %find mid-threshold value
    mid_threshold = mid_fraction * OCOG_amp;
    
    upper_bin_temp = find(bin(firstbin:lastbin) > mid_threshold) + firstbin-1;
    upper_bin = upper_bin_temp(1);
    lower_bin = upper_bin - 1;
    
    %interpolate to exact threshold
    fraction = (mid_threshold - bin(lower_bin))/(bin(upper_bin)-bin(lower_bin));
    mid_bin_position = (lower_bin-1) + fraction;
    bin_offset = win_centre - mid_bin_position;
    
    RTpoint(i) = bin_offset;
end

RTpoint = RTpoint*sampleSize;
end

    


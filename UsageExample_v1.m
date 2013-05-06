%download ExampleData.mat from my DropBox at https://www.dropbox.com/home/CRYOSAT2_FBR2SAR
%this folder also contains some pertinent papers to read FYI.
load ExampleData.mat
RETRACK_FLAG = 1;
HAMWIN_FLAG = 0;
[BURSTPOWER] = CRYOSAT2_FBR2SAR(COMPLEX_DATA_I,COMPLEX_DATA_Q,SAT_VEL,ALT_RATE,ALTITUDE,RETRACK_FLAG,HAMWIN_FLAG);
%output will be a [50x64x256] matrix of power of beams in every burst
%average across all bursts, and plot a radargram
figure;pcolor(squeeze(mean(BURSTPOWER,1)).');shading interp

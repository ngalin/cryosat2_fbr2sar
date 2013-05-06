%Interpolating and resampling function, where:
% I = interpolation rate
% D = resampling location, number of samples to shift
% x(n) = input waveform
% y(n*(1/I) + D) = output waveform

function [y] = interpResampFFT(x,winSelect,I,D,delT)

    %interpolation and resampling
    N = length(x);
    padX = [x zeros(1,N)];
    N = 2*N;
    
    fft_shiftX = fftshift(fft(padX,N));
    rotVal = exp(2*pi*sqrt(-1)*(D/delT/N).*([0:N-1]-(N/2)));
    rot_fftShiftX = fft_shiftX.*rotVal;
    y = ifft(ifftshift(rot_fftShiftX),N);
    y = y(1:N/2);
    
    %upsampling
    N = N/2;
    
    %padZeros = zeros(1,(I-1)*N-1);
    padZeros = zeros(1,(I-1)*N);
    
    fftY = fft(y,N);
    if (isempty(padZeros))
        z = [fftY(1:N)];
    else
       z = [fftY(1:N/2) padZeros fftY(N/2+1:end)];
       z(N/2+1) = z(N/2+1)/2;
       z(N+N/2+1) = z(N/2+1);
       
       %z = [fftY(1:N/2-1) fftY(N/2) padZeros fftY(N/2) fftY(N/2+1:end)];
       %z = [fftY(1:N/2) padZeros fftY(N/2+1:end)];
        
    end
    ifftZ = ifft(z);
    %rescaling
    ifftZ = ifftZ*I;
    
    y = ifftZ;
end
   

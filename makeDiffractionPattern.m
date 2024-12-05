function [amplitude, phase, diffractionPattern] = makeDiffractionPattern

    amplitude = mean(double(imread('/Users/katespada/Documents/MATLAB/GS/GS/VFSIL0069.jpg')),3);
    amplitude(amplitude<200) = 0;
    amplitude(amplitude~=0) = 255;
    amplitude = (255-amplitude)/255;

    amplitude = padarray(amplitude,[200,200],0,"both") + 0.01;

    size(amplitude)

    phase = pi*rot90(amplitude); % This is the new phase profile
    field = amplitude .* exp(1i*phase);


    diffractionPattern = abs( ifftshift( fft2( fftshift( field ) ) ) );

end
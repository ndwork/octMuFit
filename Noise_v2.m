function [noise] = Noise_v2(stdev,signalSize)
% function [noise] = snrToNoise(stdev,signalSize)
% *Creates a noise vector assuming a Rayleigh distribution in the magnitude
%   domain of FD-OCT [See William Ling's proof for why this is
%   a valid assumption.]
%
% *Assumes that variance stays constant even if signal magnitude changes
%   [See near eqn 7 in Choma "Sensitivity advantage of swept source and
%   Fourier domain optical coherence tomography"]
%
% INPUTS:
% stdev = standard deviation of the noise (inherent to the system)
%     Reasonable value for Telesto = 9.7847
%     Most clinical systems are on par with the Telesto in terms of noise.
%     A reasonable value for a low end system is 13.6831
% signalSize = size of the signal vector
%
% OUTPUTS:
% noise = noise vector with specified standard deviation

% note: scale factor for Rayleigh distribution  = sigma

  sigma = stdev / sqrt((4-pi)/2);
  noise = raylrnd(sigma,signalSize);

end

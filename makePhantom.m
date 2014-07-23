
function [I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom( ...
  phantomType, z0, noiseProportion, muValues, muThicks)

  alpha = 0.2;
  beta = 1.0;
  L0 = 2000;
  zR = 0.2;
  maxZ = 2;
  dz = 10d-3;
  z = ( 0 : dz : maxZ )';

  if nargin < 1
    phantomType = 1;
  end

  if nargin < 2
    z0 = 0.6;
  end

  if nargin < 3
    noiseProportion = 1d-5;
  end
  
  if nargin < 5

    switch phantomType
      case 0   % Uniform phantom
        muValues = 2 * ones(1,4);    % Uniform Phantom
        muThicks = [ 0.4 0.8 0.3 ];
      case 1   % Main results
        muValues = [ 0.5 1 4 2 ];
        muThicks = [ 0.4 0.8 0.3 ];
      case 2   % Bladder
        muValues = [ 0.49 2.0 1.38 ];
        muThicks = [ 0.05 0.4 ];
      case 3   % Retina
        muValues = [ 4.78 4.60 4.78 4.60 4.84 4.60 4.90 5.03 6.01 5.82 6.13 3 ];
        muThicks = [ 0.035 0.052 0.035 0.04 0.03 0.067 0.08 0.021 0.011 0.013 0.035 ];
      case 4   % Coronary Athersclerotic lesion
        muValues = [ 2 9 3 ];
        muThicks = [ 0.4 0.6 ];
    end

  end


  

  trueMu = zeros( numel(z), 1 );
  muDepths = cumsum( muThicks );
  for i=numel(muDepths):-1:1
    lowIndx = find( z >= muDepths(i), 1 );
    if i==numel(muDepths)
      trueMu(lowIndx:end) = muValues(i+1);
    else
      trueMu(lowIndx:lastIndx-1) = muValues(i+1);
    end
    lastIndx = lowIndx;
  end
  trueMu(1:lastIndx-1) = muValues(1);

  I = mu2I( trueMu, z, z0, zR, alpha, beta, L0 );

  noise = noiseProportion * L0 * randn( numel(I), 1 );
  I = I + noise;

end

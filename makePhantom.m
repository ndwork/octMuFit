
function [I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom( varargin )
  %[I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom( ...
  %  phantomType, lambda, deltaLambda, dLambda, z0, noiseProportion, ...
  %  muValues, muThicks, zR )

  defaultPhantomType = 1;

  defaultZ0 = 0.6;
  defaultNoiseProportion = 1d-5;
  defaultMuValues = [ 0.5 1 4 2 ];
  defaultMuThicks = [ 0.4 0.8 0.3 ];
  defaultZR = 0.1059 * 2 * 1.37;
  p = inputParser;
  p.addOptional( 'phantomType', defaultPhantomType, @isnumeric );
  p.addOptional( 'z0', defaultZ0, @isnumeric );
  p.addOptional( 'noiseProportion', defaultNoiseProportion, @isnumeric );
  p.addOptional( 'muValues', defaultMuValues );
  p.addOptional( 'muThicks', defaultMuThicks );
  p.addOptional( 'zR', defaultZR, @isnumeric );
  p.parse(varargin{:});
  phantomType = p.Results.phantomType;
  z0 = p.Results.z0;
  noiseProportion = p.Results.noiseProportion;
  muValues = p.Results.muValues;
  muThicks = p.Results.muThicks;
  zR = p.Results.zR;


  switch phantomType
    case 0   % Uniform phantom
      defaultMuValues = 2 * ones(1,4);    % Uniform Phantom
      defaultMuThicks = [ 0.4 0.8 0.3 ];
    case 1   % Main results
      defaultMuValues = [ 0.5 1 4 2 ];
      defaultMuThicks = [ 0.4 0.8 0.3 ];
    case 2   % Bladder
      defaultMuValues = [ 0.49 2.0 1.38 ];
      defaultMuThicks = [ 0.05 0.4 ];
    case 3   % Retina
      defaultMuValues = [ 4.78 4.60 4.78 4.60 4.84 4.60 4.90 5.03 6.01 5.82 6.13 3 ];
      defaultMuThicks = [ 0.035 0.052 0.035 0.04 0.03 0.067 0.08 0.021 0.011 0.013 0.035 ];
    case 4   % Coronary Athersclerotic lesion
      defaultMuValues = [ 2 9 3 ];
      defaultMuThicks = [ 0.4 0.6 ];
  end
  clear p;
  p = inputParser;
  p.addOptional( 'phantomType', defaultPhantomType, @isnumeric );
  p.addOptional( 'z0', defaultZ0, @isnumeric );
  p.addOptional( 'noiseProportion', defaultNoiseProportion, @isnumeric );
  p.addOptional( 'muValues', defaultMuValues );
  p.addOptional( 'muThicks', defaultMuThicks );
  p.addOptional( 'zR', defaultZR, @isnumeric );
  p.parse(varargin{:});
  z0 = p.Results.z0;
  noiseProportion = p.Results.noiseProportion;
  muValues = p.Results.muValues;
  muThicks = p.Results.muThicks;
  zR = p.Results.zR;


  alpha = 0.2;
  beta = 1.0;
  L0 = 2000;

  maxZ = 2;
  dz = 10d-3;
  z = ( 0 : dz : maxZ )';

  trueMu = muValues(1) * ones( numel(z), 1 );
  muDepths = cumsum( muThicks );
  lastIndx = numel( trueMu );
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
  I = max( I, 0 );

end

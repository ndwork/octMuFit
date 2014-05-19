
function [I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom()

  alpha = 0.2;
  beta = 1.0;
  L0 = 1.0;
  z0 = 0.6;
  zR = 0.1;
  maxZ = 2;
  dz = 10d-3;
  z = ( 0 : dz : maxZ )';
  g = 1 ./ sqrt( ((z-z0)/zR).^2 + 1 );

  %muValues = [ 0 5 1 4 2 ];
  %muThicks = [ 0.1 0.4 0.8 0.3 ];
  muValues = [ 0.5 1 4 2 ];
  muThicks = [ 0.4 0.8 0.3 ];
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

  noise = 2d-3 * randn( numel(I), 1 );
  I = I + noise;

end

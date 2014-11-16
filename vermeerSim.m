
function vermeerSim

  z0 = -1.0;
  %noiseProportion = 1d-5;
  noiseProportion = 0;
  muValues = [1 3];
  muThicks = [1];

  [I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom( ...
    1, z0, noiseProportion, muValues, muThicks);

  muVer = muFitVermeer( I, z );

  plot( z, trueMu, 'k', 'LineWidth', 2 );
  ylim([0 5]);
  hold on;
  plot( z, muVer, 'LineWidth', 2 );
  legend( 'True Mu', 'Vermeer', 'Location', 'Northwest' );

end

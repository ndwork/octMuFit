
function etbDepth = findErrorEnergyTooBigDepth( offsetThreshPercent, trueMu, muFit, z )
  diff = abs( muFit - trueMu );
  M = numel(muFit);
  errEnergy = zeros( M, 1 );
  metric = zeros(M,1);
  normTrue = norm( trueMu, 2 );
  for m=1:M
    errEnergy(m) = norm( diff(1:m), 2 );
    metric(m) = errEnergy(m) / normTrue;
    if metric(m) > offsetThreshPercent
      etbIndx = m;
      etbDepth = z( m );
      break;
    end
  end
  %metric = errEnergy ./ normTrue;
  %etbIndx = find( metric > offsetThreshPercent, 1, 'first' );
  %etbDepth = z( etbIndx );
end

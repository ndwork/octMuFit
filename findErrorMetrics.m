
function [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit,trueMu,z)
  offsetThreshPercent = 0.05;
  %offsetThreshPercent = 0.10;
  [~, N] = size(muFit);

  etbDepths = zeros(1,N);
  vDepths = zeros(1,N);
  
  for i=1:N
    thisETB = findErrorEnergyTooBigDepth( offsetThreshPercent, ...
      trueMu(:,i), muFit(:,i), z );
    thisV = findViolationDepth( offsetThreshPercent, ...
      trueMu(:,i), muFit(:,i), z );
    
    etbDepths(i) = thisETB;
    vDepths(i) = thisV;
  end
  
  meanEtbDepth = mean(etbDepths);
  meanVDepth = mean(vDepths);
end


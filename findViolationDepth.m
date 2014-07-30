
function violationDepth = findViolationDepth( offsetThreshPercent, trueMu, muFit, z )
  offsetThresh = trueMu * offsetThreshPercent;
  diff = abs( muFit - trueMu );
  %Offset by 10 to account for Gaussian Blur Kernel
  violations = diff > offsetThresh;
  violationIndx = find( violations(10:end), 1, 'first' );
  violationIndx = violationIndx + 9;
  violationDepth = z( violationIndx );
end


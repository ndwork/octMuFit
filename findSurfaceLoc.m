
function surfaceLoc = findSurfaceLoc( line )

  n = numel(line);
  sortedLine =  sort(line);
  low = mean( sortedLine(1:round(0.1*n)) );
  high = mean( sortedLine(round(0.9*n):end) );

  thresh = ( high - low ) / 2 + low;
  
  mask = line > thresh;
  locs = find( line > thresh );
  
  % Find the first point that is below the threshold and the first (after
  % the threshold) that is above the threshold  
  maskDiffs = diff(mask);
  firstAbove = 0;
  for i = 1:length(maskDiffs)
      if(maskDiffs(i) == 1)
          firstAbove = i;
          break;
      end
  end
  
  surfaceLoc = firstAbove;
  
end

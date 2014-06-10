
function [mask, noiseLevel] = makeTheMask( I, ALA )
  % Input: image matrix
  % Output: binary image of the same size as the input, where the elements
  % where mu is 0 are set to 0

  dil = 10;
  erd = 20;

  [numR, numC] = size(I);

  mask = ones(numR, numC);

  for i = 1:numC
    skinLoc = findSurfaceLoc(I(:,i));
    mask(1:skinLoc, i) = 0;
  end

  noiseLevel = median( I(mask==0) );

  if ~ALA
    mask = mask .* ( I > noiseLevel*1.5 );
  end

  mask = imdilate(mask, strel('square', dil));
  mask = imerode(mask, strel('square', erd));

end

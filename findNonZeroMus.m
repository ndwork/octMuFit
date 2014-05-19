
function mask = findNonZeroMus(image)
  % Input: image matrix
  % Output: binary image of the same size as the input, where the elements
  % where mu is 0 are set to 0
  
  dil = 10;
  erd = 20;
  
  [numR, numC] = size(image);

  mask = ones(numR, numC);

  for i = 1:numC
    skinLoc = findSkinLoc(image(:,i));
    mask(1:skinLoc, i) = 0;
  end
  
  mask = imdilate(mask, strel('square', dil));
  mask = imerode(mask, strel('square', erd));
  

end
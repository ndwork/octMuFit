
function [out1, out2, out3] = proxF1D(x1, x2, x3, mask, b, lambda, eta)

  % Note: F2 and F3 are the same function.

  out1 = proxF1(x1, mask, b, lambda);
  out2 = proxF2(x2, mask, lambda, eta);
  out3 = proxF2(x3, mask, lambda, eta);

end


function out = proxF1(x, mask, b, lambda)
  out = ( lambda ./ (1+lambda*mask) ) .* ( mask.*b + x./lambda );
end


function out = proxF2(x, mask, lambda, eta)

  out = zeros(size(x));
  
  largeXIndxs = find( mask==1 & x >= eta*lambda );
  out(largeXIndxs) = x(largeXIndxs) - eta*lambda;
  smallXIndxs = find( mask==1 & x < -eta*lambda );
  out(smallXIndxs) = x(smallXIndxs) + eta*lambda;
  
end


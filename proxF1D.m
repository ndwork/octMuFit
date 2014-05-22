function [out1, out2] = proxF1D(x1, x2, mask, b, lambda, eta)

  out1 = proxF1(x1, mask, b, lambda);
  out2 = proxF2(x2, mask, lambda, eta);

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


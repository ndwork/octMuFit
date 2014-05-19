function [out1, out2] = proxF1D(x1, x2, b, lambda, eta)

  out1 = proxF1(x1, b, lambda);
  out2 = proxF2(x2, lambda, eta);

end

function out = proxF1(x1, b, lambda)
  %out = (lambda/(lambda+1)) * (x1 + b);
  out = (x1+lambda*b)/(1+lambda);
end

function out = proxF2(x2, lambda, eta)

  out = zeros(size(x2));
  
  for i = 1:numel(x2)
    if(x2(i) >= 0)
      out(i) = max( x2(i) - lambda*eta, 0 );
    else
      out(i) = min( x2(i) + lambda*eta, 0 );
    end
  end
  
end


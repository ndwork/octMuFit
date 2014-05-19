function [out1, out2, out3] = proxF(x1, x2, x3, b, lambda, eta)

  out1 = proxF1(x1, b, lambda);
  out2 = proxF2(x2, lambda, eta);
  out3 = proxF2(x3, lambda, eta);

end

function out = proxF1(x1, b, lambda)
  out = (lambda/(lambda+1)) * (x1 + b);
end

function out = proxF2(x2, lambda, eta)

  out = zeros(size(x2));
  
  for i = 1:numel(x2)
    if(x2(i) >= 0)
      out(i) = x2(i) - lambda*eta;
    else
      out(i) = x2(i) + lambda*eta;
    end
  end
  
end
function [out1, out2] = proxFconj1D(z1, z2, b, sigma, eta)

  out1 = proxSigF1conj(z1, b, sigma);
  out2 = proxSigfF2conj(z2, eta);

end

function out = proxSigF1conj(z1, b, sigma)
  out = (z1 - sigma*b)/(1+sigma);
end

function out = proxSigfF2conj(z2, eta)
  out = min(z2, eta);
  out = max(out, -eta);
end


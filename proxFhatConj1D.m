
function [out1, out2] = proxFhatConj1D( z1, z2, bHat, sigma, eta, normA, normD )
  out1 = proxSigF1conj(z1, bHat, sigma, normA);
  out2 = proxSigfF2conj(z2, eta, normD);
end

function out = proxSigF1conj(z1, bHat, sigma, normA)
  nASq = normA*normA;
  out = 2*nASq*sigma / (sigma + 2*nASq)*(z1/sigma - bHat);
end

function out = proxSigfF2conj(z2, eta, normD)
    out = min(z2, eta*normD);
    out = max(out, -eta*normD);
end

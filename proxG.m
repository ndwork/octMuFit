function out = proxG(z, mask)
  out = max(z, 0);
  z(mask==0) = 0;
end
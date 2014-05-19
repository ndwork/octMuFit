function out = proxG(z, mask)
  out = max(z, 0);
  zeroMaskIndxs = find( mask == 0 );
  z(zeroMaskIndxs) = 0;
end
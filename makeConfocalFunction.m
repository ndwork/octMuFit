
function h = makeConfocalFunction( z, z0, zR )

  tmp = (z-z0)/zR;
  h = 1 ./ sqrt( tmp.^2 + 1 );

end
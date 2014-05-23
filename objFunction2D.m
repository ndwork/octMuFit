function out = objFunction2D(gam,I,mask,dz,dx,z,z0,zR,etaZ, etaX)
  
  b = makeb(I);
  Ag = applyA(gam, I, dz, z, z0, zR);
  Dgz = applyD(gam, dz, 1);
  Dgx = applyD(gam, dx, 2);
  
  out = 0.5 * sum ( mask(:).*(Ag(:) - b(:)).^2 ) + etaZ*sum(abs(mask(:).*Dgz(:))) + etaX*sum(abs(mask(:).*Dgx(:)));

end

% function[Einc] = PlaneWavePointWise(Eo, kvec, r)
% 
% Point-wise evaluation of plane wave at voxel centres. 
% (In contrast to integration of plane wave over voxels, as required for
% Galerkin VIE.)

function[Einc] = PlaneWavePointWise(Eo, kvec, r)

krx = kvec(1).*r(:,:,:,1);
kry = kvec(2).*r(:,:,:,2);
krz = kvec(3).*r(:,:,:,3);

kr = krx+kry+krz;

expKr=exp(1i*kr);

Einc(:,:,:,1) = Eo(1).*expKr;
Einc(:,:,:,2) = Eo(2).*expKr;
Einc(:,:,:,3) = Eo(3).*expKr;
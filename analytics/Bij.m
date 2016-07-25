function bessIJ = Bij(i,j,a,b)

Iia = besseli(i,a);
Kjb = besselk(j,b);
Ijb = besseli(j,b);
Kia = besselk(i,a);

bessIJ = Iia.*Kjb+(-1)^(i-j+1)*Ijb.*Kia;
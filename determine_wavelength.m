function wavelength = determine_wavelength(Axis,Wave)

if isrow(Axis); Axis = Axis'; end
if isrow(Wave); Wave = Wave'; end

deriv = abs(diff(Wave));
deriv = [deriv; median(deriv)];

[a,b] = sort(deriv);
[c,d] = sort(b(1:20));

x = diff(c);
y = median(x(x > 5));
wavelength = 2*(Axis(1)-Axis(y)); 

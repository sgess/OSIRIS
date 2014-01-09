function wavelength = determine_wavelength(Axis,Wave)

if isrow(Axis); Axis = Axis'; end
if isrow(Wave); Wave = Wave'; end

% chop off naughty bits
nel = numel(Wave);
maxel = round(0.9*nel);

deriv = abs(diff(Wave(1:maxel)));
deriv = [deriv; median(deriv)];

[~,b] = sort(deriv);
[c,~] = sort(b(1:20));

x = diff(c);
y = round(median(x(x > 10)));
wavelength = 2*(Axis(1)-Axis(y)); 
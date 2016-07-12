% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;
close all;
%%
cmap = custom_cmap;
savE = 0;

%% FIND DATA

data_dir = '/Users/sgess/Desktop/sims/data/os_tars/2016/';
plot_dir = '/Users/sgess/Desktop/plots/OS/';
date_dir = '2016/Jul/12/'; date_par = '2016/Jul/12/';
set_dir = 'ele131/'; plot_name = 'ele131'; 

data_loc = [data_dir set_dir];
plot_loc = [plot_dir set_dir];
paramloc = ['params/' date_par];
ext = '.eps'; ext_type = 'epsc';
if(~exist(plot_loc,'dir'))
    mkdir(plot_loc);
end

try 
    load([paramloc 'param_' set_dir(1:end-1) '.mat']);
    E0 = param_struct.plasma.field;
    skin_depth = param_struct.plasma.SD;
    n0 = param_struct.plasma.density;
catch
    disp('No param file');
    n0 = 1e17;
    [omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0] = plasma_parameters(n0);
end


%% LOAD DATA

file_number = 1;
num_str = num2str(file_number,'%06d');

bp_file = ['MS/DENSITY/beam/charge/charge-beam-' num_str '.h5'];
pp_file = ['MS/DENSITY/plasma/charge/charge-plasma-' num_str '.h5'];
e1_file = ['MS/FLD/e1/e1-' num_str '.h5'];
e2_file = ['MS/FLD/e2/e2-' num_str '.h5'];
b3_file = ['MS/FLD/b3/b3-' num_str '.h5'];

p1_file = ['MS/PHA/p1p2/beam/p1p2-beam-' num_str '.h5'];
x1_file = ['MS/PHA/p1x1/beam/p1x1-beam-' num_str '.h5'];

bp_type = 'charge';
pp_type = 'charge';
e1_type = 'e1';
e2_type = 'e2';
b3_type = 'b3';
p1_type = 'p1p2';
x1_type = 'p1x1';


% Load beam charge density
beam_rho = LOAD_DATA([data_loc bp_file],bp_type);
plas_rho = LOAD_DATA([data_loc pp_file],pp_type);
field_e1 = LOAD_DATA([data_loc e1_file],e1_type)*E0;
field_e2 = LOAD_DATA([data_loc e2_file],e2_type);
field_b3 = LOAD_DATA([data_loc b3_file],b3_type);

[z_axis, r_axis] = LOAD_AXIS([data_loc bp_file]);
zz = skin_depth*z_axis;
rr = skin_depth*r_axis;

plas_length = size(plas_rho,1);

RR = linspace(r_axis(1),r_axis(2),size(beam_rho,2));
ZZ = linspace(z_axis(1),z_axis(2),size(beam_rho,1)) - z_axis(1);

RAXIS = linspace(rr(1),rr(2),size(beam_rho,2));
ZAXIS = fliplr(linspace(zz(1),zz(2),size(beam_rho,1))-zz(1));
%% PLOT DENSITY

ions = -repmat(plas_rho(plas_length,:),plas_length,1);
charge = beam_rho'+plas_rho'+ions';

figure(1);
imagesc(ZAXIS,RAXIS,plas_rho');
axis xy;
axis image;
colormap(cmap.bwr);
caxis([-1 1]);
colorbar;
xlabel('Z [\mum]','fontsize',16);
ylabel('R [\mum]','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
title('Plasma Density','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_plasma_rho' ext],ext_type); end;


figure(2);
imagesc(ZAXIS,RAXIS,beam_rho');
axis xy;
axis image;
colormap(cmap.bwr);
caxis([-1 1]);
colorbar;
xlabel('Z [\mum]','fontsize',16);
ylabel('R [\mum]','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'n_0' ,'fontsize',16);
title('Beam Density','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_beam_rho' ext],ext_type); end;

figure(3);
imagesc(ZAXIS,RAXIS,charge);
axis xy;
axis image;
colormap(cmap.bwr);
caxis([-1 1]);
colorbar;
xlabel('Z [\mum]','fontsize',16);
ylabel('R [\mum]','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'n_0' ,'fontsize',16);
title('Charge Density','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_charge_rho' ext],ext_type); end;

%% PLOT FIELD

emax = max(abs(field_e1(:,1)));

figure(4);
imagesc(ZAXIS,RAXIS,field_e1');
axis xy;
axis image;
xlabel('Z [\mum]','fontsize',16);
ylabel('R [\mum]','fontsize',16);
colormap(cmap.bwr);
caxis([-emax emax]);
colorbar;
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'E_z (GV/m)','fontsize',16);
title('Longitudinal Field','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_EZ' ext],ext_type); end;

figure(5);
wavelength = determine_wavelength(ZAXIS,field_e1(:,1));
plot(ZAXIS,field_e1(:,1));
xlabel('Z [\mum]','fontsize',16);
ylabel('E_z (GV/m)','fontsize',16);
title(['On Axis Longitudinal Field, \lambda = ' num2str(wavelength,'%0.2f') '\mum'],'fontsize',16);
nz = length(field_e1(:,1));
dont_count = zeros(nz,1);
dont_count(1:(nz-100)) = 1;
[a,b] = max(dont_count.*field_e1(:,1));
a = double(a);
hold on;
plot(ZAXIS(b),a,'r*');
hold off;
text(ZAXIS(b),a,[' E_{max} = ' num2str(a,'%0.2f') ' GV/m']);
if savE; saveas(gca,[plot_loc plot_name '_EZ_axis' ext],ext_type); end;


% calculate fft
figure(6);
npow = nextpow2(nz);        % next power of 2
nfft = 2^npow;
y = fft(field_e1(:,1),nfft)/nz; % fft has units of 1/(# of cells)
dz = (ZAXIS(1)-ZAXIS(2))*1e-6; % intercell spacing in m
Fs = 1/dz;                     % wavenumber spacing in 1/m
f = Fs*linspace(0,1,nfft/2+1)/2; % only plot below nyquist freq (factor of 2)
realY = 2*abs(y(1:(nfft/2+1))); % plot real part of FFT
semilogy(f,realY);
[a,b] = max(realY);
hold on;
semilogy(f(b),realY(b),'r*');
hold off;
xlabel('Reduced Wavenumber (k/2 \pi) [m^{-1}]','fontsize',16);
title('Fourier transform of on axis E_z field','fontsize',16);
text_str = [' \lambda = ' num2str(1e6/f(b),'%0.2f') ' \mum'];
xpos = f(b);
ypos = double(realY(b));
text(xpos,ypos,text_str);
if savE; saveas(gca,[plot_loc plot_name '_FFT' ext],ext_type); end;

figure(7);
bfield = -(field_e2'-field_b3');
bmax = max(abs(bfield(:)));
imagesc(ZAXIS,RAXIS,bfield);
axis xy;
axis image;
xlabel('Z [\mum]','fontsize',16);
ylabel('R [\mum]','fontsize',16);
colormap(cmap.bwr);
caxis([-bmax bmax]);
colorbar;
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'E_r - B_{\theta} (MT/m)','fontsize',16);
title('Focusing Field (for Positrons)','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_ER' ext],ext_type); end;
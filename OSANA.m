% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;

savE = 1;

%data_dir = '/Users/sgess/Desktop/FACET/os_tars/';
%plot_dir = '/Users/sgess/Desktop/FACET/OS_PLOTS/';

data_dir = '/Users/sgess/Desktop/data/os_tars/';
plot_dir = '/Users/sgess/Desktop/plots/OS/';

%date_dir = '2012/Sep/07/';
%date_dir = '2013/Mar/29/';
%date_dir = '2013/May/01/';
%date_dir = '2013/May/02/';
%date_dir = '2013/May/03/';
%date_dir = '2013/May/06/';
%date_dir = '2013/Sep/04/';
date_dir = '2013/now/';

%set_dir = 'OS_eShort2/';
%set_dir  = 'hollow2/';
%set_dir  = 'test_2/';
%set_dir  = 'e_test/';
%set_dir  = 'dense_p/';
%set_dir  = 'dense_e/';
%set_dir  = 'thin_p/';
%set_dir  = 'thin_e/';
%set_dir  = 'long_test/';
%set_dir  = 'long_test2/';
%set_dir  = 'long_test3/';
%set_dir = 'long_testE/';
%set_dir = 'e_d2/';
%set_dir = 'p_r1/';
%set_dir = 'e_t1/';
%set_dir = 'p_t1/';
%set_dir = 'e_r3/';
%set_dir = 'p_r3/';
%set_dir = 'e_r3_17/';
%set_dir = 'e_r3_18/';
%set_dir = 'wtest/'; plot_name = 'sd_01';   % 0.1 sd
%set_dir = 'wtest2/'; plot_name = 'holl';  % hollow channel
%set_dir = 'wtest3/'; plot_name = 'sd_10';  % 1.0 sd
%set_dir = 'wtest4/'; plot_name = 'sd_05';  % 0.5 sd
set_dir = 'wtest5/'; plot_name = 'full';  % plasma everywhere

n0 = 1e17;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0] = plasma_parameters(n0);

data_loc = [data_dir date_dir set_dir];
plot_loc = [plot_dir date_dir set_dir];
ext = '.eps'; ext_type = 'epsc';
if(~exist(plot_loc,'dir'))
    mkdir(plot_loc);
end


for i = 1:1;

file_number = i;
num_str = num2str(file_number,'%06d');

bp_file = ['MS/DENSITY/beam/charge/charge-beam-' num_str '.h5'];
pp_file = ['MS/DENSITY/plasma/charge/charge-plasma-' num_str '.h5'];
e1_file = ['MS/FLD/e1/e1-' num_str '.h5'];
e2_file = ['MS/FLD/e2/e2-' num_str '.h5'];
b3_file = ['MS/FLD/b3/b3-' num_str '.h5'];

bp_type = 'charge';
pp_type = 'charge';
e1_type = 'e1';
e2_type = 'e2';
b3_type = 'b3';

% DNeg = [1 0 0;
%         1 1 0;
%         1 1 1;];
% 
% DPos = [1 1 1;
%         0 1 1;
%         0 0 1;];
%     
% 
D = [1 0 0; 1 1 0; 1 1 1; 0 1 1; 0 0 1;];
F = [0 0.25 0.5 0.75 1];
G = linspace(0,1,256);
cmap = interp1(F,D,G);


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

ions = -repmat(plas_rho(plas_length,:),plas_length,1);
%charge = flipdim(beam_rho'+plas_rho'+ions',1);
charge = beam_rho'+plas_rho'+ions';

zero_frac = -min(min(charge))/(max(max(charge))-min(min(charge)));

% GPos = linspace(0,1,round((1-zero_frac)*256));
% GNeg = linspace(0,0.83,round(zero_frac*256));
% 
% cpos = interp1(F,DPos,GPos);
% cneg = interp1(F,DNeg,GNeg);
% cmap = [cneg; cpos];

figure(1);
%imagesc(ZZ,RR,abs(log(flipdim(plas_rho',1))));
imagesc(ZAXIS,RAXIS,plas_rho');
axis xy;
axis image;
colormap(cmap);
caxis([-1 1]);
colorbar;
xlabel('\mum','fontsize',16);
ylabel('\mum','fontsize',16);
%xlabel('c/\omega_p','fontsize',16);
%ylabel('c/\omega_p','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', ['n_0 [' num2str(n0,'%1.1e') ']' ],'fontsize',16);
title('Plasma Density','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_plasma_rho' ext],ext_type); end;


figure(2);
%imagesc(ZZ,RR,abs(log(flipdim(beam_rho',1))));
imagesc(ZAXIS,RAXIS,beam_rho');
axis xy;
axis image;
colormap(cmap);
caxis([-1 1]);
colorbar;
xlabel('\mum','fontsize',16);
ylabel('\mum','fontsize',16);
%xlabel('c/\omega_p','fontsize',16);
%ylabel('c/\omega_p','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', ['n_0 [' num2str(n0,'%1.1e') ']' ],'fontsize',16);
title('Beam Density','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_beam_rho' ext],ext_type); end;

figure(3);
imagesc(ZAXIS,RAXIS,charge);
axis xy;
axis image;
colormap(cmap);
caxis([-1 1]);
colorbar;
xlabel('\mum','fontsize',16);
ylabel('\mum','fontsize',16);
%xlabel('c/\omega_p','fontsize',16);
%ylabel('c/\omega_p','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', ['n_0 [' num2str(n0,'%1.1e') ']' ],'fontsize',16);
title('Charge Density','fontsize',16);
if savE; saveas(gca,[plot_loc plot_name '_charge_rho' ext],ext_type); end;

% figure;
% imagesc(ZZ,RR(122)-RR(2:122),flipdim(plas_rho(:,2:122)',1));
% colorbar;
% %caxis([min(min(plas_rho(:,2:122)))/10 0]);
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
% title('Plasma Density','fontsize',16);
% v =  axis;
% text(3*v(2)/5,5*v(4)/6,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% %saveas(gca,[plot_dir date_dir set_dir 'plas_rho.pdf']);

% figure;
% imagesc(ZZ,RR(122)-RR(2:122),flipdim(beam_rho(:,2:122)',1));
% colorbar;
% %caxis([0 max(max(beam_rho(:,2:122)))/10]);
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
% title('Beam Density','fontsize',16);
% v =  axis;
% text(3*v(2)/5,5*v(4)/6,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% %saveas(gca,[plot_dir date_dir set_dir 'beam_rho.pdf']);
% 

  emax = max(abs(field_e1(:,1)));
  
  figure(4);
  imagesc(ZAXIS,RAXIS,field_e1');
  axis xy;
  axis image;
  xlabel('\mum','fontsize',16);
  ylabel('\mum','fontsize',16);
  colormap(cmap);
  caxis([-emax emax]);
  colorbar;
  t = colorbar('peer',gca);
  set(get(t,'ylabel'),'String', 'E_z (GV/m)','fontsize',16);
  title('Longitudinal Field','fontsize',16);
  if savE; saveas(gca,[plot_loc plot_name '_EZ' ext],ext_type); end;

  figure(5);
  plot(ZAXIS,field_e1(:,1));
  xlabel('\mum','fontsize',16);
  ylabel('E_z (GV/m)','fontsize',16);
  title('On Axis Longitudinal Field','fontsize',16);
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
  xlabel('\mum','fontsize',16);
  ylabel('\mum','fontsize',16);
  colormap(cmap);
  caxis([-bmax bmax]);
  colorbar;
  t = colorbar('peer',gca);
  set(get(t,'ylabel'),'String', 'E_r - B_{\theta} (MT/m)','fontsize',16);
  title('Focusing Field (for Positrons)','fontsize',16);
  if savE; saveas(gca,[plot_loc plot_name '_ER' ext],ext_type); end;

  
% imagesc(ZZ,RR(122)-RR(2:122),flipdim(field_e1(:,2:122)',1));
% colorbar;
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'm c \omega_p / e','fontsize',16);
% title('Longitudinal E Field','fontsize',16);
% v =  axis;
% text(3*v(2)/5,5*v(4)/6,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% %saveas(gca,[plot_dir date_dir set_dir 'EZ_2D.pdf']);
% 
% figure;
% plot(ZZ,field_e1(:,2));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('m c \omega_p / e','fontsize',16);
% title('Longitudinal E Field on Axis','fontsize',16);
% v =  axis;
% text(2*v(2)/3,4*v(3)/5,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% %saveas(gca,[plot_dir date_dir set_dir 'EZ_1D.pdf']);
% save('../COMPARE/OS_pShort2.mat','ZZ','field_e1');

end
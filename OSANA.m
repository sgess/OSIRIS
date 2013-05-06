% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;

savE = 0;

data_dir = '/Users/sgess/Desktop/FACET/os_tars/';
plot_dir = '/Users/sgess/Desktop/FACET/OS_PLOTS/';

%data_dir = '/Users/sgess/Desktop/data/os_tars/';
%plot_dir = '/Users/sgess/Desktop/plots/OS/';

%date_dir = '2012/Sep/07/';
%date_dir = '2013/Mar/29/';
%date_dir = '2013/May/01/';
%date_dir = '2013/May/02/';
date_dir = '2013/May/03/';

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
set_dir = 'e_t1/';


n0 = 1e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0] = plasma_parameters(n0);

data_loc = [data_dir date_dir set_dir];
plot_loc = [plot_dir date_dir set_dir];
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
charge = flipdim(beam_rho'+plas_rho'+ions',1);
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
imagesc(ZAXIS,RAXIS,flipdim(plas_rho',1));
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
if savE; saveas(gca,[plot_loc 'plasma_rho.png']); end;


figure(2);
%imagesc(ZZ,RR,abs(log(flipdim(beam_rho',1))));
imagesc(ZAXIS,RAXIS,flipdim(beam_rho',1));
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
if savE; saveas(gca,[plot_loc 'beam_rho.png']); end;

figure(3);
imagesc(ZAXIS,RAXIS,charge);
axis xy;
axis image;
colormap(cmap);
caxis([-2 2]);
colorbar;
xlabel('\mum','fontsize',16);
ylabel('\mum','fontsize',16);
%xlabel('c/\omega_p','fontsize',16);
%ylabel('c/\omega_p','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', ['n_0 [' num2str(n0,'%1.1e') ']' ],'fontsize',16);
title('Charge Density','fontsize',16);
if savE; saveas(gca,[plot_loc 'charge_rho.png']); end;

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
  figure(4);
  imagesc(ZAXIS,RAXIS,flipdim(field_e1',1));
  axis image;
  xlabel('\mum','fontsize',16);
  ylabel('\mum','fontsize',16);
  colormap(cmap);
  caxis([-0.2 0.2]);
  %caxis([-0.03 0.03]);
  colorbar;
  t = colorbar('peer',gca);
  set(get(t,'ylabel'),'String', 'E_z (GeV)','fontsize',16);
  title('Longitudinal Field','fontsize',16);
  if savE; saveas(gca,[plot_loc 'EZ.png']); end;

  figure(5);
  plot(ZAXIS,field_e1(:,1));
  xlabel('\mum','fontsize',16);
  ylabel('E_z (GeV)','fontsize',16);
  title('On Axis Longitudinal Field');
  if savE; saveas(gca,[plot_loc 'EZ_axis.png']); end;

  figure(6);
  y = fft(field_e1(:,1),512);
  Y = abs(y(1:length(y)/2));
  [max_Y,max_ind] = max(Y)
  semilogy(Y);
  
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
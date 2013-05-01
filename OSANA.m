% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;

%data_dir = '/Users/sgess/Desktop/FACET/os_tars/';
%plot_dir = '/Users/sgess/Desktop/FACET/OS_PLOTS/';

data_dir = '/Users/sgess/Desktop/data/os_tars/';
plot_dir = '/Users/sgess/Desktop/plots/OS/';

date_dir = '2013/Mar/29/';
set_dir  = 'hollow2/';

data_loc = [data_dir date_dir set_dir];
plot_loc = [plot_dir date_dir set_dir];
if(~exist(plot_loc,'dir'))
    mkdir(plot_loc);
end

for i = 0:1;

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

% Load beam charge density
beam_rho = LOAD_DATA([data_loc bp_file],bp_type);
plas_rho = LOAD_DATA([data_loc pp_file],pp_type);
field_e1 = LOAD_DATA([data_loc e1_file],e1_type);
field_e2 = LOAD_DATA([data_loc e2_file],e2_type);
field_b3 = LOAD_DATA([data_loc b3_file],b3_type);

[z_axis, r_axis] = LOAD_AXIS([data_loc bp_file]);

RR = linspace(r_axis(1),r_axis(2),length(beam_rho));
ZZ = linspace(z_axis(1),z_axis(2),length(beam_rho)) - z_axis(1);

figure(1);
%imagesc(ZZ,RR,abs(log(flipdim(plas_rho',1))));
imagesc(ZZ,RR,flipdim(plas_rho',1));
colorbar;
xlabel('c/\omega_p','fontsize',16);
ylabel('c/\omega_p','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
title('Plasma Density','fontsize',16);

figure(2);
%imagesc(ZZ,RR,abs(log(flipdim(beam_rho',1))));
imagesc(ZZ,RR,flipdim(beam_rho',1));
colorbar;
xlabel('c/\omega_p','fontsize',16);
ylabel('c/\omega_p','fontsize',16);
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
title('Beam Density','fontsize',16);

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
  figure(3);
  imagesc(ZZ,RR,flipdim(field_e1',2));
  colorbar;
  
  figure(4);
  plot(ZZ,field_e1(:,1));
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
% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;
close all;
%%
cmap = custom_cmap;
savE = 0;

%% FIND DATA

%data_dir = '/Users/sgess/Desktop/sims/data/os_tars/2016/';
%plot_dir = '/Users/sgess/Desktop/plots/OS/';

data_dir = '/Users/sgess/Desktop/sims/OS/data/2017/';
plot_dir = '/Users/sgess/Desktop/sims/OS/plots/2017/';

%date_dir = '2016/Jul/16/'; date_par = '2016/Jul/16/';
%date_dir = '2016/Aug/21/'; date_par = '2016/Aug/21/';
%date_dir = '2016/Aug/22/'; date_par = '2016/Aug/22/';
%date_dir = '2016/Aug/23/'; date_par = '2016/Aug/23/';

%date_dir = '2017/Jan/10/'; date_par = '2017/Jan/10/';
%date_dir = '2017/Jan/15/'; date_par = '2017/Jan/15/';
%date_dir = '2017/Jan/16/'; date_par = '2017/Jan/16/';
%date_dir = '2017/Jan/17/'; date_par = '2017/Jan/17/';
date_dir = '2017/Jan/22/'; date_par = '2017/Jan/22/';

%set_dir = 'ele131/'; plot_name = 'ele131';
%set_dir = 'hol131/'; plot_name = 'hol131';
%set_dir = 'cdfTest1/'; plot_name = 'cdfTest1';
%set_dir = 'cdfTest2/'; plot_name = 'cdfTest2';
%set_dir = 'moreCDF2/'; plot_name = 'moreCDF2';
%set_dir = 'width8/'; plot_name = 'width8';
%set_dir = 'width16/'; plot_name = 'width16';
%set_dir = 'width32/'; plot_name = 'width32';
%set_dir = 'yisss/'; plot_name = 'yisss';
%set_dir = 'e225_expt8/'; plot_name = 'e225_expt8';
%set_dir = 'e225_cdf1/'; plot_name = 'e225_cdf1';
%set_dir = 'e225_gauss2/'; plot_name = 'e225_gauss2';
%set_dir = 'e225_trap1/'; plot_name = 'e225_trap1';
%set_dir = 'e225_trap7/'; plot_name = 'e225_trap7';
%set_dir = 'e225_gauss3/'; plot_name = 'e225_gauss3';
%set_dir = 'e225_gauss5/'; plot_name = 'e225_gauss5';
%set_dir = 'trap_comp3/'; plot_name = 'trap_comp3';
%set_dir = 'partial1/'; plot_name = 'partial1';

%set_dir = 'trap_comp3/'; plot_name = 'trap_comp3';
%set_dir = 'partialExp/'; plot_name = 'partialExp';
%set_dir = 'partial2/'; plot_name = 'partial2';
%set_dir = 'partial3/'; plot_name = 'partial3';
%set_dir = 'partial4/'; plot_name = 'partial4';
%set_dir = 'partial5/'; plot_name = 'partial5';
%set_dir = 'partial0/'; plot_name = 'partial0';
%set_dir = 'theOnePercent/'; plot_name = 'theOnePercent';
%set_dir = 'ramp1pct1/'; plot_name = 'ramp1pct1';
%set_dir = 'ramp20pct1/'; plot_name = 'ramp20pct1';
%set_dir = 'pct1long80/'; plot_name = 'pct1long80';
%set_dir = 'pct2long80/'; plot_name = 'pct2long80';
%set_dir = 'long162part0/'; plot_name = 'long162part0';
%set_dir = 'uniFuck4/'; plot_name = 'uniFuck4';
%set_dir = 'long80part02/'; plot_name = 'long80part02';
set_dir = 'long80part02wid10/'; plot_name = 'long80part02wid10';


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
%     r_in = param_struct.hollow.radius - param_struct.hollow.width/2;
%     r_out = param_struct.hollow.radius + param_struct.hollow.width/2;
%     a = r_in/skin_depth;
%     b = r_out/skin_depth;
    N1 = param_struct.beam.N_particles;
    charge1 = param_struct.beam.charge;
    sigma_z = param_struct.beam.sigma_z;
    z_cen = param_struct.pos.beam_Z*skin_depth;
    
catch
    disp('No param file');
    n0 = 1e17;
    [omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0] = plasma_parameters(n0);
end


%% LOAD DATA

file_number = 4;
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
ZAXIS = fliplr(linspace(zz(1),zz(2),size(beam_rho,1))-(zz(2)-z_cen));
%% PLOT DENSITY

ions = -repmat(plas_rho(plas_length,:),plas_length,1);
charge = beam_rho'+plas_rho'+ions';
p_rho = plas_rho';
d_rho = p_rho(:,end);
B_rho = beam_rho';
b_rho = B_rho(1,:);

PLOT_OS2('density',ZAXIS,RAXIS,plas_rho',cmap.bwr,11);
PLOT_OS2('d1',[],RAXIS,d_rho,[],12);
PLOT_OS2('d1',[],ZAXIS,b_rho,[],13);

PLOT_OS2('density',ZAXIS,RAXIS,beam_rho',cmap.bwr,21);
PLOT_OS2('density',ZAXIS,RAXIS,charge,cmap.bwr,31,'cax',0.1);
PLOT_OS2('ez2',ZAXIS,RAXIS,field_e1',cmap.bwr,41);
PLOT_OS2('ez1',ZAXIS,[],field_e1(:,1),cmap.bwr,51);
PLOT_OS2('fr2',ZAXIS,RAXIS,-(field_e2'-field_b3'),cmap.bwr,61);

%% Compare Theory
zz = fliplr(ZAXIS);
beam_cent = ZAXIS(1) - skin_depth*param_struct.pos.beam_Z;
[EZ_out, rho_b, shift_ind] = CompareTheory(zz,field_e1(:,1)',n0,N1,charge1,sigma_z,beam_cent,a,b);

figure(10);
plot(ZAXIS,field_e1(:,1)','b',ZAXIS,EZ_out,'r');

%%
shape32 = plas_rho(end-2,:);
ez32 = field_e1(:,1);
%%
figure(1);
plot(RAXIS,-shape0,'b',RAXIS,-shape4,'c',RAXIS,-shape8,'g',RAXIS,-shape16,'m',RAXIS,-shape32,'r','linewidth',3);
axis([0 158 0 1]);
legend('\sigma = 0 \mum','\sigma = 4 \mum','\sigma = 8 \mum','\sigma = 16 \mum','\sigma = 32 \mum','location','northwest');
xlabel('R [\mum]');
ylabel('n/n_0');
set(gca,'fontsize',18);

figure(2);
plot(ZAXIS,ez0,'b',ZAXIS,ez4,'c',ZAXIS,ez8,'g',ZAXIS,ez16,'m',ZAXIS,ez32,'r','linewidth',3);
axis([0 630 -0.25 0.25]);
legend('\sigma = 0 \mum','\sigma = 4 \mum','\sigma = 8 \mum','\sigma = 16 \mum','\sigma = 32 \mum','location','southwest');
xlabel('Z [\mum]');
ylabel('E_z [GV/m]');
set(gca,'fontsize',18);
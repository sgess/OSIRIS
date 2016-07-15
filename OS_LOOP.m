% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;
close all;
cmap = custom_cmap;
savE = 0;

header = '/Users/sgess/Desktop/sims/data/os_tars/2016/';
file_number = 1;
types= {'beam_charge','plasma_charge','EZ','ER','BPHI'};

date_pars = {'2016/Jul/12/','2016/Jul/14/','2016/Jul/14/','2016/Jul/14/','2016/Jul/14/'};
set_dirs = {'hol131/','moreCDF2/','width8/','width16/','width32/'};
n_sets = length(set_dirs);

shapes = zeros(160,n_sets);
ezs = zeros(640,n_sets);

for i = 1:length(set_dirs)
    
    date_par = date_pars{i};
    set_dir = set_dirs{i};


    d_struct(i) = AG_DATA(set_dir,date_par,header,types,file_number);

    plas_rho = d_struct(i).plasma_charge';
    beam_rho = d_struct(i).beam_charge';
    ion_rho = d_struct(i).ion_charge';
    total_rho = beam_rho+plas_rho+ion_rho;
    EZ = d_struct(i).EZ';
    ER = d_struct(i).ER';
    BPHI = d_struct(i).BPHI';
    ZAXIS = d_struct(i).z_axis;
    RAXIS = d_struct(i).r_axis;
    
    EZ_LINE = EZ(1,:);
    FR = -(ER-BPHI);

    shapes(:,i) = plas_rho(:,end-2);
    ezs(:,i) = EZ_LINE;
    
    %PLOT_OS2('density',ZAXIS,RAXIS,plas_rho,cmap.bwr,1);
    %PLOT_OS2('density',ZAXIS,RAXIS,beam_rho,cmap.bwr,2);
    PLOT_OS2('density',ZAXIS,RAXIS,total_rho,cmap.bwr,3,'cax',0.1);
    PLOT_OS2('ez2',ZAXIS,RAXIS,EZ,cmap.bwr,4);
    PLOT_OS2('ez1',ZAXIS,[],EZ_LINE,cmap.bwr,5);
    %PLOT_OS2('fr2',ZAXIS,RAXIS,FR,cmap.bwr,6);
    
    pause;
end
%%
figure(1);
plot(RAXIS,-shapes(:,1),'b',RAXIS,-shapes(:,2),'c',RAXIS,-shapes(:,3),'g',RAXIS,-shapes(:,4),'m',RAXIS,-shapes(:,5),'r','linewidth',3);
axis([0 158 0 1]);
legend('\sigma = 0 \mum','\sigma = 4 \mum','\sigma = 8 \mum','\sigma = 16 \mum','\sigma = 32 \mum','location','northwest');
xlabel('R [\mum]');
ylabel('n/n_0');
set(gca,'fontsize',18);

figure(2);
plot(ZAXIS,ezs(:,1),'b',ZAXIS,ezs(:,2),'c',ZAXIS,ezs(:,3),'g',ZAXIS,ezs(:,4),'m',ZAXIS,ezs(:,5),'r','linewidth',3);
axis([0 630 -0.01 0.01]);
legend('\sigma = 0 \mum','\sigma = 4 \mum','\sigma = 8 \mum','\sigma = 16 \mum','\sigma = 32 \mum','location','southwest');
xlabel('Z [\mum]');
ylabel('E_z [GV/m]');
set(gca,'fontsize',18);

%%
beam_cent = d_struct(1).param.pos.beam_Z*d_struct(1).param.plasma.SD;
beam_z = ZAXIS(1)-beam_cent;
z_short = ZAXIS(250:end);

[a1,b1] = min(ezs(250:end,1));
z1 = z_short(b1);
d1 = z1 - beam_z;
[a2,b2] = min(ezs(250:end,2));
z2 = z_short(b2);
d2 = z2 - beam_z;
[a3,b3] = min(ezs(250:end,3));
z3 = z_short(b3);
d3 = z3 - beam_z;
[a4,b4] = min(ezs(250:end,4));
z4 = z_short(b4);
d4 = z4 - beam_z;
[a5,b5] = min(ezs(250:end,5));
z5 = z_short(b5);
d5 = z5 - beam_z;

ds = [d1 d2 d3 d4 d5];
zs = [z1 z2 z3 z4 z5];
wids = [0 4 8 16 32];

plot(wids,ds,'ko')

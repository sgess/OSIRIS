% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;
close all;
%%
cmap = custom_cmap;
savE = 0;

%% FIND DATA

data_dir = '/Users/sgess/Desktop/sims/OS/data/2017/';
plot_dir = '/Users/sgess/Desktop/sims/OS/plots/2017/';

date_dirs = {'2017/Jan/16/';
             '2017/Jan/10/';
             '2017/Jan/15/';
             '2017/Jan/15/';
             '2017/Jan/15/';
             '2017/Jan/15/'};

date_pars = {'2017/Jan/16/';
             '2017/Jan/10/';
             '2017/Jan/15/';
             '2017/Jan/15/';
             '2017/Jan/15/';
             '2017/Jan/15/'};

set_dirs = {'partial0/';
            'partial1/';
            'partial2/';
            'partial3/';
            'partial4/';
            'partial5/'};

beam_rhos = zeros(640,320,6);
plas_rhos = zeros(640,320,6);
field_e1s = zeros(640,320,6);
field_e2s = zeros(640,320,6);
field_b3s = zeros(640,320,6);

d_rhos = zeros(320,6);
ezs = zeros(640,6);

for i = 1:6
            
    set_dir = set_dirs{i};
    date_par = date_pars{i};
    
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
        r_in = param_struct.hollow.radius - param_struct.hollow.width/2;
        r_out = param_struct.hollow.radius + param_struct.hollow.width/2;
        a = r_in/skin_depth;
        b = r_out/skin_depth;
        N1 = param_struct.beam.N_particles;
        charge1 = param_struct.beam.charge;
        sigma_z = param_struct.beam.sigma_z;
        
    catch
        disp('No param file');
        n0 = 1e17;
        [omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0] = plasma_parameters(n0);
    end
    
    
    
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
    beam_rhos(:,:,i) = LOAD_DATA([data_loc bp_file],bp_type);
    plas_rhos(:,:,i) = LOAD_DATA([data_loc pp_file],pp_type);
    field_e1s(:,:,i) = LOAD_DATA([data_loc e1_file],e1_type)*E0;
    field_e2s(:,:,i) = LOAD_DATA([data_loc e2_file],e2_type);
    field_b3s(:,:,i) = LOAD_DATA([data_loc b3_file],b3_type);
    
    beam_rho = beam_rhos(:,:,i);
    plas_rho = plas_rhos(:,:,i);
    field_e1 = field_e1s(:,:,i);
    field_e2 = field_e2s(:,:,i);
    field_b3 = field_b3s(:,:,i);
    
    
    
    [z_axis, r_axis] = LOAD_AXIS([data_loc bp_file]);
    zz = skin_depth*z_axis;
    rr = skin_depth*r_axis;
    
    plas_length = size(plas_rho,1);
    
    RR = linspace(r_axis(1),r_axis(2),size(beam_rho,2));
    ZZ = linspace(z_axis(1),z_axis(2),size(beam_rho,1)) - z_axis(1);
    
    RAXIS = linspace(rr(1),rr(2),size(beam_rho,2));
    ZAXIS = fliplr(linspace(zz(1),zz(2),size(beam_rho,1))-zz(1));
    
    ions = -repmat(plas_rho(plas_length,:),plas_length,1);
    charge = beam_rho'+plas_rho'+ions';
    p_rho = plas_rho';
    d_rho = p_rho(:,end);
    B_rho = beam_rho';
    b_rho = B_rho(1,:);
    d_rhos(:,i) = d_rho;
    ez = field_e1(:,1);
    ezs(:,i) = field_e1(:,1);
    
    
%     PLOT_OS2('density',ZAXIS,RAXIS,plas_rho',cmap.bwr,11);
%     PLOT_OS2('d1',[],RAXIS,d_rho,[],12);
%     PLOT_OS2('d1',[],ZAXIS,b_rho,[],13);
%     
%     PLOT_OS2('density',ZAXIS,RAXIS,beam_rho',cmap.bwr,21);
%     PLOT_OS2('density',ZAXIS,RAXIS,charge,cmap.bwr,31,'cax',0.1);
%     PLOT_OS2('ez2',ZAXIS,RAXIS,field_e1',cmap.bwr,41);
%     PLOT_OS2('ez1',ZAXIS,[],ez,cmap.bwr,51);
%     PLOT_OS2('fr2',ZAXIS,RAXIS,-(field_e2'-field_b3'),cmap.bwr,61);
%     
%     pause;
end
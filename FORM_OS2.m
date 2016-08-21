function os_struct = FORM_OS2(param_struct)

eol = ',\n';
sep = ', ';

os_struct.simulation.n0 = [strrep(num2str(param_struct.plasma.density,'%0.1E'),'E+','d') '\n'];

os_struct.node_conf.node_number = [num2str(param_struct.sim.nodeZ) ', ' num2str(param_struct.sim.nodeR) eol];
os_struct.node_conf.if_periodic = '.false., .false.,\n';

os_struct.grid.nx_p = [num2str(param_struct.size.Grid_Z) sep num2str(param_struct.size.Grid_R) eol];

os_struct.time_step.dt = [num2str(param_struct.time.d_norm) eol];
os_struct.time_step.ndump = [num2str(param_struct.sim.ndump) eol];

os_struct.space.rmin = [num2str(0.0) sep num2str(0.0) eol];
os_struct.space.rmax = [num2str(param_struct.size.Box_Z) sep num2str(param_struct.size.Box_R) eol];
os_struct.space.if_move = '.true., .false.,\n';

os_struct.time.tmin = [num2str(0.0) sep 'tmax = ' num2str(param_struct.time.L_tot) eol];

os_struct.particles.num_species = [num2str(param_struct.sim.N_species) eol];

os_struct.species(1).name = ['''beam''' eol];
os_struct.species(1).num_par_max = [num2str(2000000) eol];
os_struct.species(1).rqm = [num2str(param_struct.beam.rqm) eol];
os_struct.species(1).num_par_x = [num2str(4) sep num2str(4) '\n'];
os_struct.species(1).vth = [num2str(0.0) sep num2str(param_struct.beam.angle_r) sep num2str(param_struct.beam.angle_r) eol];
os_struct.species(1).vfl = [num2str(0.0) sep num2str(0.0) sep num2str(0.0) eol];
os_struct.species(1).interpolation = '"quadratic"\n';
os_struct.species(1).free_stream = [param_struct.sim.free_stream eol];
os_struct.species(1).den_min = [strrep(num2str(param_struct.sim.den_min,'%0.3E'),'E','d') eol];
os_struct.species(1).dgam = [num2str(param_struct.time.d_gamma) eol];
os_struct.species(1).num_dgam = [num2str(param_struct.time.num_dgam) eol];
os_struct.species(1).num_ddgam = [num2str(param_struct.time.num_ddgam) eol];

os_struct.profile(1).profile_type = '"gaussian", "gaussian",\n';
os_struct.profile(1).gauss_center = [num2str(param_struct.pos.beam_Z) sep num2str(param_struct.pos.beam_R) eol];
os_struct.profile(1).gauss_sigma = [num2str(param_struct.beam.kpsz) sep num2str(param_struct.beam.kpsr) '\n'];
os_struct.profile(1).gauss_rangeZ = [num2str(param_struct.pos.Range_Z_min) sep num2str(param_struct.pos.Range_Z_max) eol];
os_struct.profile(1).gauss_rangeR = [num2str(0.0) sep num2str(param_struct.pos.Range_R_max) eol];
os_struct.profile(1).density = [num2str(param_struct.beam.ratio) eol];

os_struct.species(2).name = ['''plasma''' eol];
os_struct.species(2).num_par_max = [num2str(2000000) eol];
os_struct.species(2).rqm = [num2str(param_struct.plasma.rqm) eol];
os_struct.species(2).num_par_x = [num2str(4) sep num2str(4) '\n'];
os_struct.species(2).vth = [num2str(0.0) sep num2str(0.0) eol];
os_struct.species(2).vfl = [num2str(0.0) sep num2str(0.0) eol];
os_struct.species(2).interpolation = '"quadratic"\n';

if param_struct.hollow.use == 0
    os_struct.profile(2).num_x = [num2str(6) eol];
    % plasma N z
    os_struct.profile(2).fx1 = [num2str(0.0) sep num2str(0.0) sep num2str(1.0)...
        sep num2str(1.0) sep num2str(1.0) sep num2str(1.0) eol];
    % plasma Z z
    os_struct.profile(2).x1 = [num2str(0.0) sep num2str(param_struct.pos.plasma_Z_start)...
        sep num2str(param_struct.pos.plasma_Z_ramp) sep num2str(param_struct.pos.plasma_Z_end/2)...
        sep num2str(3*param_struct.pos.plasma_Z_end/4) sep num2str(param_struct.pos.plasma_Z_end) eol];
    % plasma N r
    os_struct.profile(2).fx2 = [num2str(1.0) sep num2str(1.0) sep num2str(1.0)...
        sep num2str(1.0) sep num2str(0.0) sep num2str(0.0) eol];
    % plasma R r
    os_struct.profile(2).x2 = [num2str(0.0) sep num2str(param_struct.size.Box_R/2)...
        sep num2str(3*param_struct.size.Box_R/4) sep num2str(param_struct.pos.plasma_R_end)...
        sep num2str(param_struct.pos.plasma_R_ramp) sep num2str(param_struct.size.Box_R) eol];
else
    if strcmp(param_struct.hollow.type,'flat') && param_struct.hollow.n_points == 6
        os_struct.profile(2).num_x = [num2str(param_struct.hollow.n_points) eol];
        % plasma N z
        os_struct.profile(2).fx1 = [num2str(0.0) sep num2str(0.0) sep num2str(1.0)...
            sep num2str(1.0) sep num2str(1.0) sep num2str(1.0) eol];
        % plasma Z z
        os_struct.profile(2).x1 = [num2str(0.0) sep num2str(param_struct.pos.plasma_Z_start)...
            sep num2str(param_struct.pos.plasma_Z_ramp) sep num2str(param_struct.pos.plasma_Z_end/2)...
            sep num2str(3*param_struct.pos.plasma_Z_end/4) sep num2str(param_struct.pos.plasma_Z_end) eol];
        % plasma N r
        os_struct.profile(2).fx2 = [num2str(param_struct.hollow.N_vec(1)) sep num2str(param_struct.hollow.N_vec(2))...
            sep num2str(param_struct.hollow.N_vec(3)) sep num2str(param_struct.hollow.N_vec(4))...
            sep num2str(param_struct.hollow.N_vec(5)) sep num2str(param_struct.hollow.N_vec(6)) eol];
        % plasma R r
        os_struct.profile(2).x2 = [num2str(param_struct.hollow.R_vec(1)) sep num2str(param_struct.hollow.R_vec(2))...
            sep num2str(param_struct.hollow.R_vec(3)) sep num2str(param_struct.hollow.R_vec(4))...
            sep num2str(param_struct.hollow.R_vec(5)) sep num2str(param_struct.hollow.R_vec(6)) eol];
    elseif strcmp(param_struct.hollow.type,'hollow') && param_struct.hollow.n_points == 6
        os_struct.profile(2).num_x = [num2str(param_struct.hollow.n_points) eol];
        % plasma N z
        os_struct.profile(2).fx1 = [num2str(0.0) sep num2str(0.0) sep num2str(1.0)...
            sep num2str(1.0) sep num2str(1.0) sep num2str(1.0) eol];
        % plasma Z z
        os_struct.profile(2).x1 = [num2str(0.0) sep num2str(param_struct.pos.plasma_Z_start)...
            sep num2str(param_struct.pos.plasma_Z_ramp) sep num2str(param_struct.pos.plasma_Z_end/2)...
            sep num2str(3*param_struct.pos.plasma_Z_end/4) sep num2str(param_struct.pos.plasma_Z_end) eol];
        % plasma N r
        os_struct.profile(2).fx2 = [num2str(param_struct.hollow.N_vec(1)) sep num2str(param_struct.hollow.N_vec(2))...
            sep num2str(param_struct.hollow.N_vec(3)) sep num2str(param_struct.hollow.N_vec(4))...
            sep num2str(param_struct.hollow.N_vec(5)) sep num2str(param_struct.hollow.N_vec(6)) eol];
        % plasma R r
        os_struct.profile(2).x2 = [num2str(param_struct.hollow.R_vec(1)) sep num2str(param_struct.hollow.R_vec(2))...
            sep num2str(param_struct.hollow.R_vec(3)) sep num2str(param_struct.hollow.R_vec(4))...
            sep num2str(param_struct.hollow.R_vec(5)) sep num2str(param_struct.hollow.R_vec(6)) eol];
    elseif strcmp(param_struct.hollow.type,'cdf')
        os_struct.profile(2).num_x = [num2str(param_struct.hollow.n_points) eol];
        os_struct.profile(2).nx = num2str(param_struct.hollow.n_points);
        
        % N z string 
        n_z_start = [num2str(0.0,'%0.0f') sep num2str(0.0,'%0.0f')];
        n_z_mid = [];
        for i = 1:(param_struct.hollow.n_points-2); n_z_mid = [n_z_mid sep num2str(1.0,'%0.0f')]; end;
        n_z_end = eol;
        os_struct.profile(2).fx1 = [n_z_start n_z_mid n_z_end];
        
        % Z z string
        z_vals = linspace(param_struct.pos.plasma_Z_ramp,param_struct.pos.plasma_Z_end,param_struct.hollow.n_points);
        z_z_start = [num2str(0.0,'%0.0f') sep num2str(param_struct.pos.plasma_Z_start,'%0.0f') sep num2str(param_struct.pos.plasma_Z_ramp,'%0.0f')];
        z_z_mid = [];
        for i = 4:numel(z_vals); z_z_mid = [z_z_mid sep num2str(z_vals(i),'%0.0f')]; end;
        z_z_end = eol;
        os_struct.profile(2).x1 = [z_z_start z_z_mid z_z_end];
        
        % N r string
        n_r_start = num2str(param_struct.hollow.N_vec(1),'%0.2f');
        n_r_mid = [];
        for i = 2:param_struct.hollow.n_points; n_r_mid = [n_r_mid sep num2str(param_struct.hollow.N_vec(i),'%0.2f')]; end;
        n_r_end = eol;
        os_struct.profile(2).fx2 = [n_r_start n_r_mid n_r_end];
        
        % R r string
        r_r_start = num2str(param_struct.hollow.R_vec(1),'%0.2f');
        r_r_mid = [];
        for i = 2:param_struct.hollow.n_points; r_r_mid = [r_r_mid sep num2str(param_struct.hollow.R_vec(i),'%0.2f')]; end;
        r_r_end = eol;
        os_struct.profile(2).x2 = [r_r_start r_r_mid r_r_end];
        
    elseif strcmp(param_struct.hollow.type,'cdf2')
        os_struct.profile(2).num_x = [num2str(param_struct.hollow.n_points) eol];
        os_struct.profile(2).nx = num2str(param_struct.hollow.n_points);
        
        % N z string 
        n_z_start = [num2str(0.0,'%0.0f') sep num2str(0.0,'%0.0f')];
        n_z_mid = [];
        for i = 1:(param_struct.hollow.n_points-2); n_z_mid = [n_z_mid sep num2str(1.0,'%0.0f')]; end;
        n_z_end = eol;
        os_struct.profile(2).fx1 = [n_z_start n_z_mid n_z_end];
        
        % Z z string
        z_vals = linspace(param_struct.pos.plasma_Z_ramp,param_struct.pos.plasma_Z_end,param_struct.hollow.n_points);
        z_z_start = [num2str(0.0,'%0.0f') sep num2str(param_struct.pos.plasma_Z_start,'%0.0f') sep num2str(param_struct.pos.plasma_Z_ramp,'%0.0f')];
        z_z_mid = [];
        for i = 4:numel(z_vals); z_z_mid = [z_z_mid sep num2str(z_vals(i),'%0.0f')]; end;
        z_z_end = eol;
        os_struct.profile(2).x1 = [z_z_start z_z_mid z_z_end];
        
        % N r string
        n_r_start = num2str(param_struct.hollow.N_vec(1),'%0.2f');
        n_r_mid = [];
        for i = 2:param_struct.hollow.n_points; n_r_mid = [n_r_mid sep num2str(param_struct.hollow.N_vec(i),'%0.2f')]; end;
        n_r_end = eol;
        os_struct.profile(2).fx2 = [n_r_start n_r_mid n_r_end];
        
        % R r string
        r_r_start = num2str(param_struct.hollow.R_vec(1),'%0.2f');
        r_r_mid = [];
        for i = 2:param_struct.hollow.n_points; r_r_mid = [r_r_mid sep num2str(param_struct.hollow.R_vec(i),'%0.2f')]; end;
        r_r_end = eol;
        os_struct.profile(2).x2 = [r_r_start r_r_mid r_r_end];
    else
        error('need more code');
    end
end


function os_struct = FORM_OS(param_struct)

eol = ',\n';
sep = ', ';

os_struct.simulation.n0 = [strrep(num2str(param_struct.plasma.density,'%0.3E'),'E+','d') '\n'];

os_struct.node_conf.node_number = [num2str(param_struct.sim.nodeZ) ', ' num2str(param_struct.sim.nodeX) eol];
os_struct.node_conf.if_periodic = '.false., .false.,\n';

os_struct.grid.nx_p = [num2str(param_struct.size.Grid_Z) sep num2str(param_struct.size.Grid_X) eol];
os_struct.grid.coordinates = [param_struct.sim.coordinates eol];

os_struct.time_step.dt = [num2str(param_struct.time.d_norm) eol];
os_struct.time_step.ndump = [num2str(param_struct.sim.ndump) eol];

os_struct.space.xmin = [num2str(0.0) sep num2str(0.0) eol];
os_struct.space.xmax = [num2str(param_struct.size.Box_Z) sep num2str(param_struct.size.Box_X) eol];
os_struct.space.if_move = '.true., .false.,\n';

os_struct.time.tmin = [num2str(0.0) sep 'tmax = ' num2str(param_struct.time.total_steps) eol];

os_struct.particles.num_species = [num2str(param_struct.sim.N_species) eol];

os_struct.species(1).name = ['''beam''' eol];
os_struct.species(1).num_par_max = [num2str(2000000) eol];
os_struct.species(1).rqm = [num2str(param_struct.beam.rqm) eol];
os_struct.species(1).num_par_x = [num2str(4) sep num2str(4) '\n'];
os_struct.species(1).vth = [num2str(0.0) sep num2str(param_struct.beam.angle_x) sep num2str(param_struct.beam.angle_y) eol];
os_struct.species(1).vfl = [num2str(0.0) sep num2str(0.0) sep num2str(0.0) eol];
os_struct.species(1).interpolation = '"quadratic"\n';
os_struct.species(1).free_stream = [param_struct.sim.free_stream eol];
os_struct.species(1).den_min = [strrep(num2str(param_struct.sim.den_min,'%0.3E'),'E','d') eol];
os_struct.species(1).dgam = [num2str(param_struct.time.d_gamma) eol];
os_struct.species(1).num_dgam = [num2str(param_struct.time.num_dgam) eol];
os_struct.species(1).num_ddgam = [num2str(param_struct.time.num_dgam-param_struct.time.d_gam_norm) eol];

os_struct.profile(1).profile_type = '"gaussian", "gaussian",\n';
os_struct.profile(1).gauss_center = [num2str(param_struct.pos.beam_Z_norm) sep num2str(param_struct.pos.beam_X_norm) eol];
os_struct.profile(1).gauss_sigma = [num2str(param_struct.beam.kpsz) sep num2str(param_struct.beam.kpsx) '\n'];
os_struct.profile(1).gauss_rangeZ = [num2str(param_struct.pos.beam_Z_norm-param_struct.pos.beam_Z_range)...
    sep num2str(param_struct.pos.beam_Z_norm+param_struct.pos.beam_Z_range) eol];
os_struct.profile(1).gauss_rangeR = [num2str(0.0) sep num2str(param_struct.pos.beam_X_norm+param_struct.pos.beam_X_range) eol];
os_struct.profile(1).density = [num2str(param_struct.beam.ratio) eol];

os_struct.species(2).name = ['''plasma''' eol];
os_struct.species(2).num_par_max = [num2str(2000000) eol];
os_struct.species(2).rqm = [num2str(param_struct.plasma.rqm) eol];
os_struct.species(2).num_par_x = [num2str(4) sep num2str(4) '\n'];
os_struct.species(2).vth = [num2str(0.0) sep num2str(0.0) eol];
os_struct.species(2).vfl = [num2str(0.0) sep num2str(0.0) eol];
os_struct.species(2).interpolation = '"quadratic"\n';

os_struct.profile(2).num_x = [num2str(6) eol];
os_struct.profile(2).fx1 = [num2str(0.0) sep num2str(0.0) sep num2str(1.0)...
    sep num2str(1.0) sep num2str(1.0) sep num2str(1.0) eol];
os_struct.profile(2).x1 = [num2str(0.0) sep num2str(param_struct.pos.plasma_Z_start)...
    sep num2str(param_struct.pos.plasma_Z_ramp) sep num2str(param_struct.pos.plasma_Z_end/2)...
    sep num2str(3*param_struct.pos.plasma_Z_end/4) sep num2str(param_struct.pos.plasma_Z_end) eol];
os_struct.profile(2).fx2 = [num2str(1.0) sep num2str(1.0) sep num2str(1.0)...
    sep num2str(1.0) sep num2str(0.0) sep num2str(0.0) eol];
os_struct.profile(2).x2 = [num2str(0.0) sep num2str(param_struct.size.Box_X/2)...
    sep num2str(3*param_struct.size.Box_X/4) sep num2str(param_struct.pos.plasma_X_end)...
    sep num2str(param_struct.pos.plasma_X_ramp) sep num2str(param_struct.size.Box_X) eol];







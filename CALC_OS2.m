function param_struct = CALC_OS2(input_struct)

% import standard SI constants
SI_consts;


%%%%%%%%%%%%%%%%%%%%%
% PLASMA ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%

param_struct.plasma.density  = input_struct.plasma.density; % plasma density [cm^-3]
param_struct.plasma.charge   = input_struct.plasma.charge;  % plasma particle charge [e]
param_struct.plasma.mass     = input_struct.plasma.mass;    % plasma particle mass [e_m]

% Calc plasma parameters
param_struct.plasma.omega_p  = sqrt(param_struct.plasma.density*...
    1e6*SI_e^2/(SI_em*SI_eps0));                                                    % plasma frequency [rad/s]
param_struct.plasma.lambda_p = 2*pi*SI_c*1e6/param_struct.plasma.omega_p;           % plasma wavelength [um]
param_struct.plasma.k_p      = param_struct.plasma.omega_p/SI_c;                    % plasma wavenumber [1/m]
param_struct.plasma.SD       = 1e6*SI_c/param_struct.plasma.omega_p;                % plasma skin depth [um]
param_struct.plasma.time     = 1/param_struct.plasma.omega_p;                       % characteristic time scale [s]
param_struct.plasma.rqm      = param_struct.plasma.mass/param_struct.plasma.charge; % mass to charge ratio
param_struct.plasma.field    = SI_em*SI_c*param_struct.plasma.omega_p/(1e9*SI_e);   % GV/m


%%%%%%%%%%%%%%%%%%%
% BEAM ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

param_struct.beam.charge      = input_struct.beam.charge;                        % beam particle charge [e]
param_struct.beam.mass        = input_struct.beam.mass;                          % beam particle mass [e_m]
param_struct.beam.N_particles = input_struct.beam.N_particles;                   % number of beam particles
param_struct.beam.rqm         = param_struct.beam.mass/param_struct.beam.charge; % mass to charge ratio
param_struct.beam.gamma       = input_struct.beam.gamma;
param_struct.beam.energy      = param_struct.beam.gamma*param_struct.beam.mass*SI_eM/1e3;

% Betatron frequecy, wavelength
param_struct.beam.omega_b  = sqrt(param_struct.plasma.density*1e6*SI_e^2/... 
    (2*param_struct.beam.gamma*param_struct.beam.mass*SI_em*SI_eps0));       % betatron frequency [rad/s]
param_struct.beam.lambda_b = 2*pi*SI_c*1e6/param_struct.beam.omega_b;        % betatron wavelength [um]

% Beam size
param_struct.beam.sigma_r = input_struct.beam.sigma_r;
param_struct.beam.sigma_z = input_struct.beam.sigma_z; 

% Radial emittance
param_struct.beam.emit_r = input_struct.beam.emit_r;

% Normalized beam divergence
param_struct.beam.angle_r = param_struct.beam.emit_r/param_struct.beam.sigma_r; % emit_x/sigma_x

% Calc beam density
param_struct.beam.density = (10000^3)*param_struct.beam.N_particles/...
    ((2*pi)^(3/2)*param_struct.beam.sigma_r^2*param_struct.beam.sigma_z);        % beam density [cm^-3]

% Calc N_b/N_p
param_struct.beam.ratio = param_struct.beam.density/param_struct.plasma.density; % beam density

% Calc k_p*sigma_x,z
param_struct.beam.kpsr = param_struct.beam.sigma_r/param_struct.plasma.SD; % beam length R [skin depths]
param_struct.beam.kpsz = param_struct.beam.sigma_z/param_struct.plasma.SD; % beam length Z [skin depths]

% Calc peak gaussian current
param_struct.beam.I_peak = 1e3*param_struct.beam.N_particles*...
    SI_e*SI_c/(sqrt(2*pi)*param_struct.beam.sigma_z);            % beam current [kA]

% Calc max bubble radius
param_struct.plasma.R_bubble = (1/0.84)*2*param_struct.plasma.SD*...
    sqrt(param_struct.beam.I_peak/17);                               % plasma bubble radius[um]


%%%%%%%%%%%%%%%%%%%
% SIZE ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

param_struct.size.cell_frac = input_struct.size.cell;                             % Cell size [skin depths]
param_struct.size.cell_size = param_struct.size.cell_frac*param_struct.plasma.SD; % Cell size [um]

% Determine simulation box size
param_struct.size.Box_R = input_struct.size.Box_R;
param_struct.size.Box_Z = input_struct.size.Box_Z;

% Determine simulation box size
param_struct.size.R = param_struct.size.Box_R*param_struct.plasma.SD; % Box size X [um]
param_struct.size.Z = param_struct.size.Box_Z*param_struct.plasma.SD; % Box size Z [um]

% Determine number of cells,
param_struct.size.Grid_R = ceil(param_struct.size.Box_R/param_struct.size.cell_frac); % Number of grid points
param_struct.size.Grid_Z = ceil(param_struct.size.Box_Z/param_struct.size.cell_frac); % Number of grid points


%%%%%%%%%%%%%%%%%%%
% TIME ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

% Test Courant Condition
CC = 0.6/sqrt(2*param_struct.size.cell_frac^-2);
if input_struct.sim.dt > CC
    error(['dt does not satisfy Courant Condition. dt must be less than' num2str(CC)]);
else
    param_struct.time.d_norm = input_struct.sim.dt;                          % Time step [omega_p^-1]
    param_struct.time.d_real = input_struct.sim.dt*param_struct.plasma.time; % Time step [s]
end

% Calculate length of simulation in skin depth units
param_struct.time.L_real  = input_struct.sim.prop*1e4;                               % Propagation length in plasma [um]
param_struct.time.L_norm  = param_struct.time.L_real/param_struct.plasma.SD;         % Propagation length in plasma [skin depths]
%param_struct.time.L_time  = ceil(param_struct.time.L_norm/param_struct.time.d_norm); % Propagation time in plasma [1/omega_p]

% Calculate total simulation time including beam initialization
param_struct.time.d_gam_norm  = input_struct.sim.gamma_steps;                           % Number of steps to accelerate beam
param_struct.time.d_gamma     = param_struct.beam.gamma/param_struct.time.d_gam_norm;   % Acceleration in gamma per step
param_struct.time.num_dgam    = ceil(param_struct.size.Box_Z/param_struct.time.d_norm); % Number of time steps before plasma
param_struct.time.num_ddgam   = param_struct.time.num_dgam - param_struct.time.d_gam_norm; % Number of non-accelerating steps
param_struct.time.L_tot       = param_struct.size.Box_Z+param_struct.time.L_norm;    % Total time of simulation [1/omega_p]



%%%%%%%%%%%%%%%%%%%%%%%
% POSITION ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%

% Calculate where to put beam
param_struct.pos.beam_Z = input_struct.pos.Center_Z;
param_struct.pos.beam_R = input_struct.pos.Center_R;
param_struct.pos.Range_R_max = input_struct.pos.Range_R_max;
param_struct.pos.Range_R_min = input_struct.pos.Range_R_min;
param_struct.pos.Range_Z_max = input_struct.pos.Range_Z_max;
param_struct.pos.Range_Z_min = input_struct.pos.Range_Z_min;

% Calculate where to put plasma
param_struct.pos.plasma_Z_start = param_struct.time.num_dgam*param_struct.time.d_norm+param_struct.size.Box_Z; % Plasma start [skin depths]
param_struct.pos.plasma_Z_ramp  = input_struct.sim.plasma_Z_ramp+param_struct.pos.plasma_Z_start;              % Plasma ramp [skin depths]
param_struct.pos.plasma_Z_end   = 10*param_struct.time.L_tot;                                                  % Plasma end [skin depths]

param_struct.hollow.use = input_struct.hollow.use;
if param_struct.hollow.use
    param_struct.hollow.n_points = input_struct.hollow.n_points;
    param_struct.hollow.type     = input_struct.hollow.type;
    param_struct.hollow.radius   = input_struct.hollow.radius;
    param_struct.hollow.width    = input_struct.hollow.width;
    param_struct.hollow.ramp     = input_struct.hollow.ramp;
    param_struct.hollow.r_ramp   = input_struct.hollow.r_ramp;
    [param_struct.hollow.R_vec, param_struct.hollow.N_vec] = hollow_channel(param_struct);
else
    param_struct.pos.plasma_R_start = 0.0;                                                                         % Plasma start [skin depths]
    param_struct.pos.plasma_R_end   = param_struct.size.Box_R-0.1;                                                 % Plasma end [skin depths]
    param_struct.pos.plasma_R_ramp  = param_struct.pos.plasma_R_end+input_struct.sim.plasma_R_ramp;                % Plasma ramp [skin depths]
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate number of nodes
if rem(param_struct.size.Grid_R,input_struct.sim.nodeR) ~= 0
    error('Grid X not divisable by transverse nodes');
elseif param_struct.size.Grid_R/input_struct.sim.nodeR < 20
    error('Too few X grids per tranverse node. Must be at least 20');
else
    param_struct.sim.nodeR = input_struct.sim.nodeR;
end

if rem(param_struct.size.Grid_Z,input_struct.sim.nodeZ) ~= 0
    error('Grid Z not divisable by longitudinal nodes');
elseif param_struct.size.Grid_Z/input_struct.sim.nodeZ < 20
    error('Too few Z grids per longitudinal node. Must be at least 20');
else
    param_struct.sim.nodeZ = input_struct.sim.nodeZ;
end

param_struct.sim.tot_steps   = ceil(param_struct.time.L_tot/param_struct.time.d_norm);
param_struct.sim.ndump       = floor(param_struct.sim.tot_steps/input_struct.sim.n_dumps);
param_struct.sim.N_species   = input_struct.sim.N_species;
param_struct.sim.free_stream = input_struct.sim.free_stream;
param_struct.sim.den_min     = input_struct.sim.density_res;


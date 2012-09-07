function param_struct = CALC_OS(input_struct)

% import standard SI constants
eval(['run ' pwd '/SI_consts.m']);



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



%%%%%%%%%%%%%%%%%%%
% BEAM ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

param_struct.beam.charge      = input_struct.beam.charge;                        % beam particle charge [e]
param_struct.beam.mass        = input_struct.beam.mass;                          % beam particle mass [e_m]
param_struct.beam.N_particles = input_struct.beam.N_particles;                   % number of beam particles
param_struct.beam.rqm         = param_struct.beam.mass/param_struct.beam.charge; % mass to charge ratio

% Calc gamma if not specified
if input_struct.beam.gamma == 0
   if input_struct.beam.energy == 0
      error('Must specify either gamma or energy');
   end
   param_struct.beam.gamma  = input_struct.beam.energy/...
       (input_struct.beam.mass*SI_eM/1e3);                  % beam gamma
   param_struct.beam.energy = input_struct.beam.energy;     % beam energy [GeV]
end

% Calc energy if not specified
if input_struct.beam.energy == 0
   if input_struct.beam.gamma == 0
      error('Must specify either gamma or energy');
   end
   param_struct.beam.energy = input_struct.beam.gamma*...
       input_struct.beam.mass*SI_eM/1e3;                    % beam energy [GeV]
   param_struct.beam.gamma  = input_struct.beam.gamma;      % beam gamma
end

% Betatron frequecy, wavelength
param_struct.beam.omega_b  = sqrt(param_struct.plasma.density*1e6*SI_e^2/... 
    (2*param_struct.beam.gamma*param_struct.plasma.mass*SI_em*SI_eps0));     % betatron frequency [rad/s]
param_struct.beam.lambda_b = 2*pi*SI_c*1e6/param_struct.beam.omega_b;        % betatron wavelength [um]

% Calc beam size if beam_match == 1
if( input_struct.beam.beam_match )
  param_struct.beam.sigma_x = 1e3*sqrt(input_struct.beam.emit_x/...
      param_struct.plamsa.k_p*sqrt(2/param_struct.beam.gamma));      % beam size X [um]
  
  param_struct.beam.sigma_y = 1e3*sqrt(input_struct.beam.emit_y/...
      param_struct.plasma.k_p*sqrt(2/param_struct.beam.gamma));      % beam size Y [um]
  
  param_struct.beam.emit_x  = input_struct.beam.emit_x;              % beam norm X emmitance [mm*mrad]
  param_struct.beam.emit_y  = input_struct.beam.emit_x;              % beam norm Y emmitance [mm*mrad]

% Calc emittance if emit_match == 1
elseif( input_struct.beam.emit_match )
  param_struct.beam.emit_x  = 1e-6*param_struct.plasma.k_p*...
      (input_struct.beam.sigma_x)^2*sqrt(param_struct.beam.gamma/2); % beam norm X emmitance [mm*mrad]
  param_struct.beam.emit_y  = 1e-6*param_struct.plasma.k_p*...
      (input_struct.beam.sigma_y)^2*sqrt(param_struct.beam.gamma/2); % beam norm Y emmitance [mm*mrad]
  
  param_struct.beam.sigma_x = input_struct.beam.sigma_x;             % beam size X [um]
  param_struct.beam.sigma_y = input_struct.beam.sigma_y;             % beam size Y [um]

else
  param_struct.beam.sigma_x = input_struct.beam.sigma_x;             % beam size X [um]
  param_struct.beam.sigma_y = input_struct.beam.sigma_y;             % beam size Y [um]
  
  param_struct.beam.emit_x  = input_struct.beam.emit_x;              % beam norm X emmitance [mm*mrad]
  param_struct.beam.emit_y  = input_struct.beam.emit_x;              % beam norm Y emmitance [mm*mrad]
end

% Normalized beam divergence
param_struct.beam.angle_x = param_struct.beam.emit_x/param_struct.beam.sigma_x; % emit_x/sigma_x
param_struct.beam.angle_y = param_struct.beam.emit_y/param_struct.beam.sigma_y; % emit_x/sigma_x

% Calc sigma_z if z_match == 1
if( input_struct.beam.z_match )
   param_struct.beam.sigma_z = 1e6*sqrt(2)*param_struct.plasma.SD;   % beam size Z [um]
else
   param_struct.beam.sigma_z = input_struct.beam.sigma_z;            % beam size Z [um]
end

% Calc beam density
param_struct.beam.density = (10000^3)*param_struct.beam.N_particles/...
    ((2*pi)^(3/2)*param_struct.beam.sigma_x*param_struct.beam.sigma_y*...
    param_struct.beam.sigma_z);                                           % beam density [cm^-3]

% Calc N_b/N_p
param_struct.beam.ratio = param_struct.beam.density/param_struct.plasma.density; % beam density

% Calc k_p*sigma_x,y,z
param_struct.beam.kpsx = param_struct.beam.sigma_x/param_struct.plasma.SD; % beam length X [skin depths]
param_struct.beam.kpsy = param_struct.beam.sigma_y/param_struct.plasma.SD; % beam length Y [skin depths]
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
param_struct.size.Box_X = ceil(input_struct.size.X_beamradii*param_struct.beam.kpsx); % Box size X [skin depths]
param_struct.size.Box_Z = ceil(2*pi*input_struct.size.Z_wavelengths);                 % Box size Z [skin depths]

% Determine simulation box size
param_struct.size.X = param_struct.size.Box_X*param_struct.plasma.SD; % Box size X [um]
param_struct.size.Z = param_struct.size.Box_Z*param_struct.plasma.SD; % Box size Z [um]

% Determine number of cells,
param_struct.size.Grid_X = ceil(param_struct.size.Box_X/param_struct.size.cell_frac); % Number of grid points
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
param_struct.time.L_steps = ceil(param_struct.time.L_norm/param_struct.time.d_norm); % Time steps in plasma

% Calculate total simulation time including beam initialization
param_struct.time.d_gam_norm  = input_struct.sim.gamma_steps;                           % Number of steps to accelerate beam
param_struct.time.d_gamma     = param_struct.beam.gamma/param_struct.time.d_gam_norm;   % Acceleration in gamma per step
param_struct.time.num_dgam    = ceil(param_struct.size.Box_Z/param_struct.time.d_norm); % Number of steps before plasma
param_struct.time.total_steps = param_struct.time.num_dgam+param_struct.time.L_steps;   % Total number of time steps in sim



%%%%%%%%%%%%%%%%%%%%%%%
% POSITION ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%

% Calculate where to put beam
param_struct.pos.beam_Z_norm = input_struct.size.Z_center*param_struct.size.Box_Z;  % Beam position from left of box [skin depths]
param_struct.pos.beam_Z_real = param_struct.pos.beam_Z_norm*param_struct.plasma.SD; % Beam position from left of box [um]
param_struct.pos.beam_X_norm = input_struct.size.X_center*param_struct.size.Box_X;  % Beam position from bottom of box [skin depths]
param_struct.pos.beam_X_real = param_struct.pos.beam_X_norm*param_struct.plasma.SD; % Beam position from bottom of box [um]

% Calculate Gaussian range
if param_struct.pos.beam_X_norm+input_struct.beam.N_sigma_x*...
        param_struct.beam.kpsx > param_struct.size.Box_X-0.1
    error('Beam larger than box! Try a smaller N_sigma value');
else
    param_struct.pos.beam_X_range = input_struct.beam.N_sigma_x*param_struct.beam.kpsx; % Beam range [skin depths]
end

if param_struct.pos.beam_Z_norm+input_struct.beam.N_sigma_z*param_struct.beam.kpsz > param_struct.size.Box_Z-0.1 ||...
        param_struct.pos.beam_Z_norm-input_struct.beam.N_sigma_z*param_struct.beam.kpsz < 0.1        
    error('Beam larger than box! Try a smaller N_sigma value');
else
    param_struct.pos.beam_Z_range = input_struct.beam.N_sigma_z*param_struct.beam.kpsz; % Beam range [skin depths]
end

% Calculate where to put plasma
param_struct.pos.plasma_Z_start = param_struct.time.num_dgam*param_struct.time.d_norm+param_struct.size.Box_Z; % Plasma start [skin depths]
param_struct.pos.plasma_Z_ramp  = input_struct.sim.plasma_Z_ramp+param_struct.pos.plasma_Z_start;              % Plasma ramp [skin depths]
param_struct.pos.plasma_Z_end   = 10*param_struct.time.total_steps;                                            % Plasma end [skin depths]
param_struct.pos.plasma_X_start = 0.0;                                                                         % Plasma start [skin depths]
param_struct.pos.plasma_X_end   = param_struct.size.Box_Z-0.1;                                                 % Plasma end [skin depths]
param_struct.pos.plasma_X_ramp  = param_struct.pos.plasma_X_end+input_struct.sim.plasma_X_ramp;                % Plasma ramp [skin depths]



%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate number of nodes
if rem(param_struct.size.Grid_X,input_struct.sim.nodeX) ~= 0
    error('Grid X not divisable by transverse nodes');
elseif param_struct.size.Grid_X/input_struct.sim.nodeX < 20
    error('Too few X grids per tranverse node. Must be at least 20');
else
    param_struct.sim.nodeX = input_struct.sim.nodeX;
end

if rem(param_struct.size.Grid_Z,input_struct.sim.nodeZ) ~= 0
    error('Grid Z not divisable by longitudinal nodes');
elseif param_struct.size.Grid_Z/input_struct.sim.nodeZ < 20
    error('Too few Z grids per longitudinal node. Must be at least 20');
else
    param_struct.sim.nodeZ = input_struct.sim.nodeZ;
end

param_struct.sim.coordinates = input_struct.sim.coordinates;
param_struct.sim.ndump       = floor(param_struct.time.total_steps/input_struct.sim.n_dumps);
param_struct.sim.N_species   = input_struct.sim.N_species;
param_struct.sim.free_stream = input_struct.sim.free_stream;
param_struct.sim.den_min     = input_struct.sim.density_res;


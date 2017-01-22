 % OSIRIS Matlab os-stdin generation example script
% S. Gessner Sep 6, 2012

clear all;

% import standard SI constants
SI_consts;

% specify template
osinput_template_file = [pwd '/osinputs/os-stdin_template'];

% specify output
timezone = 'America/Los_Angeles';
date_dir = GET_DATE_DIR(timezone);

osinput_dir = [pwd '/osinputs/' date_dir];
if ~exist(osinput_dir,'dir')
    mkdir(osinput_dir);
end

param_dir = [pwd '/params/' date_dir];
if ~exist(param_dir,'dir')
    mkdir(param_dir);
end

command_dir = [pwd '/commands/' date_dir];
if ~exist(command_dir,'dir')
    mkdir(command_dir);
end

osinput_output_name = 'uniFuck5';
osinput_output_file = [osinput_dir 'os-stdin_' osinput_output_name];

write = 1;
% check to see if you want to overwrite file
if exist(osinput_output_file,'file')
   reply = input(['File ' osinput_output_file ' exists. \n Do you want to overwrite? y/n '], 's');
   if (strcmp(reply,'n'))
      disp('Ok. That''s cool.');
      write = 0;
   end
end


% INPUT TO OSINPUT (ASSUMING NO IONS, GAUSSIAN BEAMS, AND CYLINDRICAL COORDINATES)

% simulation parameters
input_struct.sim.N_species     = 2;               % Number of particle species
input_struct.sim.dt            = 0.016;           % Time step in 1/omega_p, must satisfy courant condition
input_struct.sim.prop          = 1.000;           % propagation length of the beam [cm]
input_struct.sim.gamma_steps   = 400;             % number of time steps for beam to accelerate during initialization
input_struct.sim.plasma_R_ramp = 0.02;            % transverse plasma ramp in skin depths (to avoid noise at boundary) 
input_struct.sim.plasma_Z_ramp = 3;               % longitudinal plasma ramp in skin depths (to avoid trapped charge) 
input_struct.sim.density_res   = 1e-8;            % resolution of beam and plasma density relative to n0 
input_struct.sim.nodeR         = 2;               % number of "transverse nodes"
input_struct.sim.nodeZ         = 4;               % number of "longitudinal nodes" 
input_struct.sim.n_dumps       = 5;               % number of data dumps per run 
input_struct.sim.free_stream   = '.true.';        % false : beam evolves, true : free stream, no evolution

% plasma parameters
input_struct.plasma.density    = 1.4e15;          % /cm^3
input_struct.plasma.charge     = -1.0;            % e 
input_struct.plasma.mass       = SI_eM/SI_eM;     % Particle mass in units of electron mass

% hollow channel profile
input_struct.hollow.use        = 0;               % use hollow channel?
input_struct.hollow.radius     = 250;             % central radius in microns
input_struct.hollow.width      = 30;              % annulus width in microns
input_struct.hollow.ramp       = 0.1;             % ramp length in microns
input_struct.hollow.r_ramp     = 10;               % ramp length in microns
input_struct.hollow.n_points   = 6;               % number of points in profile
input_struct.hollow.n_in       = 0.02;            % density of partial channel
input_struct.hollow.type       = 'partial';       % profile type

% beam parameters
input_struct.beam.charge       = 1.0;             % -1 for electron, +1 for positron
input_struct.beam.mass         = SI_eM/SI_eM;     % Particle mass in units of electron mass
input_struct.beam.N_particles  = 2.50e9;          % Number of beam particles
input_struct.beam.gamma        = 40000;           % relativistic factor gamma, if 0 energy specified below
input_struct.beam.sigma_r      = 20.;             % Gaussian sigma_x [um]
input_struct.beam.sigma_z      = 320.;            % Gaussian sigma_z [um]
input_struct.beam.emit_r       = 5.0;             % normalized X emittance [mm*mrad i.e. 1e-6]

% grid size parameters
input_struct.size.cell         = 0.05;            % cell size as a fraction of the skin depth
input_struct.size.Box_R        = 8;              % box size in skin depths
input_struct.size.Box_Z        = 16;              % box size in skin depths

% grid position parameters
input_struct.pos.Center_R        = 0;             % beam centroid position in skin depths (0 is on axis)
input_struct.pos.Center_Z        = 10;            % beam centroid position in
input_struct.pos.Range_R_max     = 3.25;          % Max extent of beam in R
input_struct.pos.Range_R_min     = 0;             % Min extent of beam in R
input_struct.pos.Range_Z_max     = 15.5;          % Max extent of beam in Z
input_struct.pos.Range_Z_min     = 4.5;           % Min extent of beam in Z

% PARAMETER CALCULATOR
param_struct = CALC_OS2(input_struct);

% OSSTDIN FORMATER
os_struct = FORM_OS2(param_struct);

% OSSTDIN WRITER
if write
    WRITE_OS2(osinput_template_file, osinput_output_file, os_struct);
    save([param_dir 'param_' osinput_output_name '.mat'], 'param_struct');
    %run_dir = WRITE_CMD(command_dir, osinput_output_name, param_struct.comp.mem,...
    %    param_struct.comp.tasks, param_struct.comp.run_time);
end

%exit;

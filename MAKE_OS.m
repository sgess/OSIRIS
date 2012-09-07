% OSIRIS Matlab os-stdin generation example script
% S. Gessner Sep 6, 2012

clear all;

% import standard SI constants
eval(['run ' pwd '/SI_consts.m']);

% specify template
osinput_template_file = [pwd '/osinputs/os-stdin_template'];

% specify output
datestr = date;
[day, monthstr]  = strtok(datestr,'-');
[month, yearstr] = strtok(monthstr,'-');
[year, eonstr]   = strtok(yearstr,'-');

osinput_dir = [pwd '/osinputs/' year '/' month '/' day '/'];
if ~exist(osinput_dir,'dir')
    mkdir(osinput_dir);
end

param_dir = [pwd '/params/' year '/' month '/' day '/'];
if ~exist(param_dir,'dir')
    mkdir(param_dir);
end

% CHANGE ME!!!!!!!!!!!!!!
osinput_output_name = 'eShort2';

osinput_output_file = [osinput_dir 'os-stdin_' osinput_output_name];

write = 1;
% check to see if you want to overwrite file
if exist(osinput_output_file,'file')
   reply = input(['File ' osinput_output_file ' exists. \n Do you want to overwrite? y/n '], 's');
   if (strcmp(reply,'n'))
      disp('Ok. Forget it.');
      write = 0;
   end
end


% INPUT TO OSINPUT (ASSUMING NO IONS, GAUSSIAN BEAMS, AND CYLINDRICAL COORDINATES)

% simulation parameters
input_struct.sim.coordinates     = '"cylindrical"'; % Z,R is order of coordinates
input_struct.sim.N_species       = 2;               % Number of particle species
input_struct.sim.dt              = 0.016;           % Time step in 1/omega_p, must satisfy courant condition
input_struct.sim.prop            = 0.040;           % propagation length of the beam [cm]
input_struct.sim.gamma_steps     = 400;             % number of time steps for beam to accelerate during initialization
input_struct.sim.plasma_X_ramp   = 0.05;            % transverse plasma ramp in skin depths (to avoid noise at boundary) 
input_struct.sim.plasma_Z_ramp   = 3;               % longitudinal plasma ramp in skin depths (to avoid trapped charge) 
input_struct.sim.density_res     = 1e-8;            % resolution of beam and plasma density relative to n0 
input_struct.sim.nodeX           = 8;               % number of "transverse nodes"
input_struct.sim.nodeZ           = 8;               % number of "longitudinal nodes" 
input_struct.sim.n_dumps         = 50;              % number of data dumps per run 
input_struct.sim.free_stream     = '.false.';       % false : beam evolves, true : free stream, no evolution

% plasma parameters
input_struct.plasma.density      = 5e16;        % /cm^3
input_struct.plasma.charge       = -1.0;        % -1 for electron, +1 for positron
input_struct.plasma.mass         = SI_eM/SI_eM; % Particle mass in units of electron mass

% beam parameters
input_struct.beam.charge         = -1.0;        % -1 for electron, +1 for positron
input_struct.beam.mass           = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.beam.N_particles    = 2.0e10;      % Number of beam particles
input_struct.beam.gamma          = 40000;       % relativistic factor gamma, if 0 energy specified below
input_struct.beam.energy         = 0;           % beam mean energy [GeV], if 0 use gamma to calculate energy
input_struct.beam.sigma_x        = 60;          % Gaussian sigma_x [um]
input_struct.beam.sigma_y        = 60;          % Gaussian sigma_y [um]
input_struct.beam.sigma_z        = 14;          % Gaussian sigma_z [um]
input_struct.beam.N_sigma_x      = 5;           % N sigma range of distribtuion in X
input_struct.beam.N_sigma_y      = 5;           % N sigma range of distribtuion in Y
input_struct.beam.N_sigma_z      = 5;           % N sigma range of distribtuion in Z
input_struct.beam.emit_x         = 50.0;        % normalized X emittance [mm*mrad i.e. 1e-6]
input_struct.beam.emit_y         = 50.0;        % normalized Y emittance [mm*mrad i.e. 1e-6]
input_struct.beam.beam_match     = 0;           % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
input_struct.beam.emit_match     = 0;           % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
input_struct.beam.z_match        = 0;           % 1: override sigma_z with sqrt(2)/k_p, 0: do nothing

% size parameters
input_struct.size.cell           = 0.05;        % cell size as a fraction of the skin depth
input_struct.size.Z_wavelengths  = 2.5;         % set box length by number of plasma wavelengths
input_struct.size.X_beamradii    = 6;           % set box length by number of beam radii
input_struct.size.Z_center       = 0.8;         % Z center of beam measured from end of box
input_struct.size.X_center       = 0;           % R center of beam measured from longitudinal axis

% PARAMETER CALCULATOR
param_struct = CALC_OS(input_struct);

% OSSTDIN FORMATER
os_struct = FORM_OS(param_struct);

% OSSTDIN WRITER
if write
    WRITE_OS(osinput_template_file, osinput_output_file, os_struct);
    save([param_dir 'param_' osinput_output_name '.mat'], 'param_struct');
end
exit;

% OSIRIS ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;

data_dir = '/Users/sgess/Desktop/FACET/os_tars/2012/Sep/07/OS_pShort2/';

beam_rho_file = 'MS/DENSITY/beam/charge/charge-beam-003128.h5';
plas_rho_file = 'MS/DENSITY/plasma/charge/charge-plasma-003128.h5';
field_e1_file = 'MS/FLD/e1/e1-003128.h5';
field_e2_file = 'MS/FLD/e2/e2-003128.h5';
field_b3_file = 'MS/FLD/b3/b3-003128.h5';

% Load beam charge density
beam_rho = LOAD_DATA([data_dir beam_rho_file]);
plas_rho = LOAD_DATA([data_dir plas_rho_file]);
field_e1 = LOAD_DATA([data_dir field_e1_file]);
field_e2 = LOAD_DATA([data_dir field_e2_file]);
field_b3 = LOAD_DATA([data_dir field_b3_file]);

[z_axis, r_axis] = LOAD_AXIS([data_dir beam_rho_file]);

figure;
imagesc(abs(log(flipdim(plas_rho',1))));
colorbar;

figure;
imagesc(abs(log(flipdim(beam_rho',1))));
colorbar;

figure;
imagesc(flipdim(plas_rho',1));
colorbar;

figure;
imagesc(flipdim(beam_rho',1));
colorbar;

figure;
imagesc(flipdim(field_e1',1));
colorbar;

figure;
plot(field_e1(:,1));


function d_struct = AG_DATA(set_dir,date_par,header,types,file_number)

% Data and param dirs
data_loc = [header set_dir];
paramloc = ['params/' date_par];
num_str = num2str(file_number,'%06d');
n_type = numel(types);

% Load Param
load([paramloc 'param_' set_dir(1:end-1) '.mat']);
d_struct.param = param_struct;

% Load Axis
file = ['MS/DENSITY/beam/charge/charge-beam-' num_str '.h5'];
[z_axis, r_axis] = LOAD_AXIS([data_loc file]);
data = LOAD_DATA([data_loc file],'charge');
skin_depth = d_struct.param.plasma.SD;
zz = skin_depth*z_axis;
rr = skin_depth*r_axis;
RAXIS = linspace(rr(1),rr(2),size(data,2));
ZAXIS = fliplr(linspace(zz(1),zz(2),size(data,1))-zz(1));
d_struct.r_axis = RAXIS;
d_struct.z_axis = ZAXIS;

% Load Rest of data
for i = 1:n_type
    
    if strcmp(types{i},'beam_charge')
        type = 'charge';
        file = ['MS/DENSITY/beam/charge/charge-beam-' num_str '.h5'];
        data = LOAD_DATA([data_loc file],type);
        d_struct.beam_charge = data;
    elseif strcmp(types{i},'plasma_charge')
        type = 'charge';
        file = ['MS/DENSITY/plasma/charge/charge-plasma-' num_str '.h5'];
        data = LOAD_DATA([data_loc file],type);
        d_struct.plasma_charge = data;
        
        plas_length = size(data,1);
        ions = -repmat(data(plas_length,:),plas_length,1);
        d_struct.ion_charge = ions;
    elseif strcmp(types{i},'EZ')
        type = 'e1';
        file = ['MS/FLD/e1/e1-' num_str '.h5'];
        data = LOAD_DATA([data_loc file],type);
        d_struct.EZ = data;
    elseif strcmp(types{i},'ER')
        type = 'e2';
        file = ['MS/FLD/e2/e2-' num_str '.h5'];
        data = LOAD_DATA([data_loc file],type);
        d_struct.ER = data;
    elseif strcmp(types{i},'BPHI')
        type = 'b3';
        file = ['MS/FLD/b3/b3-' num_str '.h5'];
        data = LOAD_DATA([data_loc file],type);
        d_struct.BPHI = data;
    end
end
        







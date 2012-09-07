function [z_axis, r_axis] = LOAD_AXIS(file)

z_axis = h5read(file,'/AXIS/AXIS1'); % longitudinal axis start and stop in skin depths
r_axis = h5read(file,'/AXIS/AXIS2'); % transverse axis start and stop in skin depths
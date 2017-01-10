function [R_vec, N_vec] = hollow_channel(param_struct)

n_points = param_struct.hollow.n_points;
type = param_struct.hollow.type;

R_vec = zeros(1,n_points);
N_vec = zeros(1,n_points);

radius = param_struct.hollow.radius/param_struct.plasma.SD;
width = param_struct.hollow.width/param_struct.plasma.SD;
ramp = param_struct.hollow.ramp/param_struct.plasma.SD;
r_max = param_struct.size.Box_R;

if strcmp(type,'flat')
    
    if n_points ~= 6
        error('More code needed');
    end
    
    R_vec(1) = 0.0;
    N_vec(1) = 0.0;
    R_vec(2) = radius - width/2 - ramp;
    N_vec(2) = 0.0;
    R_vec(3) = radius - width/2;
    N_vec(3) = 1.0;
    R_vec(4) = radius + width/2;
    N_vec(4) = 1.0;
    R_vec(5) = radius + width/2 + ramp;
    N_vec(5) = 0.0;
    R_vec(6) = param_struct.size.Box_R;
    N_vec(6) = 0.0;
    
end

if strcmp(type,'hollow')
    
    if n_points ~= 6
        error('More code needed');
    end
    
    R_vec(1) = 0.0;
    N_vec(1) = 0.0;
    R_vec(2) = radius - ramp;
    N_vec(2) = 0.0;
    R_vec(3) = radius;
    N_vec(3) = 1.0;
    R_vec(4) = radius + width/2;
    N_vec(4) = 1.0;
    R_vec(5) = radius + width/2 + ramp;
    N_vec(5) = 1.0;
    R_vec(6) = param_struct.size.Box_R;
    N_vec(6) = 1.0;
    
end

if strcmp(type,'cdf')
    
    R_vec = linspace(0,param_struct.size.Box_R,n_points);
    dR_vec = R_vec(2)-R_vec(1);
    gauss = (sqrt(2*pi)*width)^(-1)*exp(-(R_vec-radius).^2/(2*width^2));
    N_vals = 1000*dR_vec*cumsum(gauss);
    N_int_vals = floor(N_vals);
    N_vec = N_int_vals/1000;
    %N_vec(end) = 0;
    
end

if strcmp(type,'cdf2')
    w = width;
    wm = w*(n_points-3);
    rs = 0:w:wm;
    m = mean(rs);
    rg = rs+(radius-m);
    R_vec = [0 rg r_max];
    
    dR_vec = [0 diff(R_vec)];
    gauss = (sqrt(2*pi)*width)^(-1)*exp(-(R_vec-radius).^2/(2*width^2));
    N_vals = 1000*cumsum(dR_vec.*gauss);
    N_int_vals = floor(N_vals);
    N_vec = N_int_vals/1000;
    %N_vec(end) = 0;
    
end

if strcmp(type,'gauss')
    
    if n_points ~= 13
        error('More code needed');
    end
    
    R_vec = [0 (radius-3.5*width) (radius-2.5*width) (radius-1.5*width) (radius-0.75*width) (radius-0.35*width) radius...
        (radius+0.35*width) (radius+0.75*width) (radius+1.5*width) (radius+2.5*width) (radius+3.5*width) r_max];
    
    %dR_vec = [0 diff(R_vec)];
    gauss = exp(-(R_vec-radius).^2/(2*width^2));
    N_vals = 1000*gauss;
    N_int_vals = floor(N_vals);
    N_vec = N_int_vals/1000;
    %N_vec(end) = 0;
    
end

if strcmp(type,'trap')
    
    if n_points ~= 6
        error('More code needed');
    end
    r_ramp = param_struct.hollow.r_ramp/param_struct.plasma.SD;
    
    r_ll = radius-width/2-r_ramp;
    r_lh = radius-width/2;
    r_hh = radius+width/2;
    r_hl = radius+width/2+r_ramp;
    
    R_vec = [0 r_ll r_lh r_hh r_hl r_max];
    N_vec = [0 0 1 1 0 0];
  
    
end

if strcmp(type,'partial')
    
    n_in = param_struct.hollow.n_in;
    r_ramp = param_struct.hollow.r_ramp/param_struct.plasma.SD;
    
    r_ll = radius-width/2-r_ramp;
    r_lh = radius-width/2;
    r_hh = radius+width/2;
    r_hl = radius+width/2+r_ramp;
    
    R_vec = [0 r_ll r_lh r_hh r_hl r_max];
    N_vec = [n_in n_in 1 1 0 0];
  
    
end
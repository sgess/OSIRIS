function WRITE_OS2(myfilein, myfileout, os_struct)
% Generate rpinput file by replacing modified lines in a standard os-stdin
% A bunch of nastiness, Don't touch!

fid = fopen(myfilein, 'r');
fidout = fopen(myfileout, 'w');

section_simulation   = 0;
section_node_conf    = 0;
section_grid         = 0;
section_time_step    = 0;
section_space        = 0;
section_time         = 0;
section_particles    = 0;
section_species      = 0;
section_profile      = 0;

n_species = 1;

row = 0;
n_row = 0;

while( row ~= -1 )

  n_row = n_row + 1;
  row = fgets(fid);
  
  if( row ~= -1 )
    % TAG SECTIONS
    if( ~isempty(strfind(row, 'simulation')) )
      section_simulation = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_simulation = 0;
    end% if
    if( ~isempty(strfind(row, 'node_conf')) )
      section_node_conf = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_node_conf = 0;
    end% if
    if( ~isempty(strfind(row, 'grid')) )
      section_grid = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_grid = 0;
    end% if
    if( ~isempty(strfind(row, 'time_step')) )
      section_time_step = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_time_step = 0;
    end% if
    if( ~isempty(strfind(row, 'space')) )
      section_space = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_space = 0;
    end% if
    if( ~isempty(strfind(row, 'time')) )
      section_time = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_time = 0;
    end% if
    if( ~isempty(strfind(row, 'particles')) )
      section_particles = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_particles = 0;
    end% if
    if( ~isempty(strfind(row, 'species')) )
      section_species = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
      section_species = 0;
    end% if
    if( ~isempty(strfind(row, 'profile')) )
      section_profile = 1;
    end% if
    if( ~isempty(strfind(row, '}')) )
        if section_profile == 1
            n_species = n_species + 1;
        end
        section_profile = 0;    
    end% if
    

    % SIMULATION SECTION
    if(~isempty(strfind(row, 'n0')) && section_simulation == 1)
       fprintf(fidout, ['  n0 = ' os_struct.simulation.n0]);

    % NODE CONF SECTION
    elseif(~isempty(strfind(row, 'node_number')) && section_node_conf == 1)
       fprintf(fidout, ['  node_number(1:2) = ' os_struct.node_conf.node_number]);
    elseif(~isempty(strfind(row, 'if_periodic') ) && section_node_conf == 1)
       fprintf(fidout, ['  if_periodic(1:2) = ' os_struct.node_conf.if_periodic]);

    % GRID SECTION
    elseif(~isempty(strfind(row, 'nx_p')) && section_grid == 1)
       fprintf(fidout, ['  nx_p(1:2)   = ' os_struct.grid.nx_p]);
       
    % TIME STEP SECTION
    elseif(~isempty(strfind(row, 'dt')) && section_time_step == 1)
       fprintf(fidout, ['  dt    = ' os_struct.time_step.dt]);
    elseif(~isempty(strfind(row, 'ndump')) && section_time_step == 1)
       fprintf(fidout, ['  ndump = ' os_struct.time_step.ndump]);
       
    % SPACE SECTION
    elseif(~isempty(strfind(row, 'xmin')) && section_space == 1)
       fprintf(fidout, ['  xmin(1:2) = ' os_struct.space.rmin]);
    elseif(~isempty(strfind(row, 'xmax')) && section_space == 1)
       fprintf(fidout, ['  xmax(1:2) = ' os_struct.space.rmax]);
    elseif(~isempty(strfind(row, 'if_move')) && section_space == 1)
       fprintf(fidout, ['  if_move   = ' os_struct.space.if_move]);
       
    % TIME SECTION
    elseif(~isempty(strfind(row, 'tmin') ) && section_time == 1)
       fprintf(fidout, ['  tmin = ' os_struct.time.tmin]);   
       
    % PARTICLES SECTION
    elseif(~isempty(strfind(row, 'num_species')) && section_particles == 1)
       fprintf(fidout, ['  num_species = ' os_struct.particles.num_species]);
 
    % SPECIES SECTION
    elseif(~isempty(strfind(row, 'name')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  name           = ' os_struct.species(1).name]);
        elseif n_species == 2
            fprintf(fidout, ['  name           = ' os_struct.species(2).name]);
        end
    elseif(~isempty(strfind(row, 'num_par_max')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  num_par_max    = ' os_struct.species(1).num_par_max]);
        elseif n_species == 2
            fprintf(fidout, ['  num_par_max    = ' os_struct.species(2).num_par_max]);
        end
    elseif(~isempty(strfind(row, 'rqm')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  rqm            = ' os_struct.species(1).rqm]);
        elseif n_species == 2
            fprintf(fidout, ['  rqm            = ' os_struct.species(2).rqm]);
        end
    elseif(~isempty(strfind(row, 'num_par_x')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  num_par_x(1:2) = ' os_struct.species(1).num_par_x]);
        elseif n_species == 2
            fprintf(fidout, ['  num_par_x(1:2) = ' os_struct.species(2).num_par_x]);
        end  
    elseif(~isempty(strfind(row, 'vth')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  vth(1:3)       = ' os_struct.species(1).vth]);
        elseif n_species == 2
            fprintf(fidout, ['  vth(1:2)       = ' os_struct.species(2).vth]);
        end
    elseif(~isempty(strfind(row, 'vfl')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  vfl(1:3)       = ' os_struct.species(1).vfl]);
        elseif n_species == 2
            fprintf(fidout, ['  vfl(1:2)       = ' os_struct.species(2).vfl]);
        end
    elseif(~isempty(strfind(row, 'interpolation')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  interpolation  = ' os_struct.species(1).interpolation]);
        elseif n_species == 2
            fprintf(fidout, ['  interpolation  = ' os_struct.species(2).interpolation]);
        end
    elseif(~isempty(strfind(row, 'free_stream')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  free_stream    = ' os_struct.species(1).free_stream]);
        end
    elseif(~isempty(strfind(row, 'den_min')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  den_min        = ' os_struct.species(1).den_min]);
        end
    elseif(~isempty(strfind(row, ' dgam')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  dgam           = ' os_struct.species(1).dgam]);
        end
    elseif(~isempty(strfind(row, 'num_dgam')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  num_dgam       = ' os_struct.species(1).num_dgam]);
        end
    elseif(~isempty(strfind(row, 'num_ddgam')) && section_species == 1)
        if n_species == 1
            fprintf(fidout, ['  num_ddgam      = ' os_struct.species(1).num_ddgam]);
        end
        
    % PROFILE SECTION
    elseif(~isempty(strfind(row, 'profile_type')) && section_profile == 1)
        if n_species == 1
            fprintf(fidout, ['  profile_type(1:2)  = ' os_struct.profile(1).profile_type]);
        end
    elseif(~isempty(strfind(row, 'gauss_center')) && section_profile == 1)
        if n_species == 1
            fprintf(fidout, ['  gauss_center(1:2)  = ' os_struct.profile(1).gauss_center]);
        end
    elseif(~isempty(strfind(row, 'gauss_sigma')) && section_profile == 1)
        if n_species == 1
            fprintf(fidout, ['  gauss_sigma(1:2)   = ' os_struct.profile(1).gauss_sigma]);
        end
    elseif(~isempty(strfind(row, 'gauss_range(1:2,1)')) && section_profile == 1)
        if n_species == 1
            fprintf(fidout, ['  gauss_range(1:2,1) = ' os_struct.profile(1).gauss_rangeZ]);
        end
    elseif(~isempty(strfind(row, 'gauss_range(1:2,2)')) && section_profile == 1)
        if n_species == 1
            fprintf(fidout, ['  gauss_range(1:2,2) = ' os_struct.profile(1).gauss_rangeR]);
        end
    elseif(~isempty(strfind(row, 'density')) && section_profile == 1)
        if n_species == 1
            fprintf(fidout, ['  density            = ' os_struct.profile(1).density]);
        end
    elseif(~isempty(strfind(row, 'num_x')) && section_profile == 1)
        if n_species == 2
            fprintf(fidout, ['  num_x     = ' os_struct.profile(2).num_x]);
        end
    elseif(~isempty(strfind(row, 'fx(1:6,1)')) && section_profile == 1)
        if n_species == 2
            fprintf(fidout, ['  fx(1:' os_struct.profile(2).nx ',1) = ' os_struct.profile(2).fx1]);
        end
    elseif(~isempty(strfind(row, 'x(1:6,1)')) && section_profile == 1)
        if n_species == 2
            fprintf(fidout, ['   x(1:' os_struct.profile(2).nx ',1) = ' os_struct.profile(2).x1]);
        end
    elseif(~isempty(strfind(row, 'fx(1:6,2)')) && section_profile == 1)
        if n_species == 2
            fprintf(fidout, ['  fx(1:' os_struct.profile(2).nx ',2) = ' os_struct.profile(2).fx2]);
        end
    elseif(~isempty(strfind(row, 'x(1:6,2)')) && section_profile == 1)
        if n_species == 2
            fprintf(fidout, ['   x(1:' os_struct.profile(2).nx ',2) = ' os_struct.profile(2).x2]);
        end
    else
       fprintf(fidout, row);
    end%if
  end%if
end%while
fclose(fid);
fclose(fidout);
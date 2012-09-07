function data = LOAD_DATA(file)

if ~isempty(strfind(file,'charge-'))
    loc = '/charge';
elseif ~isempty(strfind(file,'e1-'))
    loc = '/e1';
elseif ~isempty(strfind(file,'e2-'))
    loc = '/e2';
elseif ~isempty(strfind(file,'e3-'))
    loc = '/e3';
elseif ~isempty(strfind(file,'b1-'))
    loc = '/b1';
elseif ~isempty(strfind(file,'b2-'))
    loc = '/b2';
elseif ~isempty(strfind(file,'b3-'))
    loc = '/b3';
end

data = h5read(file,loc); 
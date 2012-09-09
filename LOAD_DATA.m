function data = LOAD_DATA(file,loc)

data = h5read(file,['/' loc]); 
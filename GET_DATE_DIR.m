function date_dir = GET_DATE_DIR()

datstr = date;
[day, monstr] = strtok(datstr,'-');
[mon, yerstr] = strtok(monstr,'-');
[yer, ~]      = strtok(yerstr,'-');
date_dir      = [yer '/' mon '/' day '/'];
function date_dir = GET_DATE_DIR(time_zone)


if nargin == 1
    dt = datetime('now','TimeZone',time_zone,'Format','dd-MMM-yyy');
    datstr = char(dt);
else
    datstr = date;
end

[day, monstr] = strtok(datstr,'-');
[mon, yerstr] = strtok(monstr,'-');
[yer, ~]      = strtok(yerstr,'-');
date_dir      = [yer '/' mon '/' day '/'];
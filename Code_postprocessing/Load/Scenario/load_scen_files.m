function files = load_scen_files(dirs)

dir_main = dirs.dir_main;
dir_scen = dirs.dir_scen;
dir_case = dirs.dir_case;

files.plot_parameters      = [dir_main,dir_scen,dir_case,'/parameter_files/plot_parameters.m'];
files.slip_rate_parameters = [dirs.dir_raid,dirs.dir_scen,dirs.dir_case, 'parameter_files/slip_rate_parameters.m'];
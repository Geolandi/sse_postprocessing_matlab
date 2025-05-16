% clear workspace
clear variables
% close all figures
close all

% add code path
addpath(genpath('./Code_postprocessing/'));

% initialize directories where data and results are
dirs.dir_main = [pwd,'/'];
dirs.dir_raid = dirs.dir_main;
dirs.dir_data = 'Website/Data/Cascadia/';
dirs.dir_scen = 'Website/Scenarios/Cascadia/';
dirs.dir_case  = 'SSE/ES/';

% select day to process
day2process = '2025-01-26';
input_dir = [dirs.dir_raid,dirs.dir_scen,dirs.dir_case, ...
    'results/',day2process,'/'];

mm2m = 1e-3;
km2m = 1e3;

% assign names to parameter files needed later (plot_parameters and
% slip_rate_parameters)
files = load_scen_files(dirs);

% load data used as input for the vbICA (i.e., after linear trend and
% offset corrections)
data_X = load([input_dir, 'X.mat']);
Xd = data_X.Xd;
clear data_X;
% load ICA results
data_ICA = load([input_dir, 'ICA.mat']);
ICA = data_ICA.ICA_essential;
ICA.type     = Xd.type;
ICA.llh      = Xd.llh;
ICA.timeline = Xd.timeline;
ICA.decmode  = Xd.decmode;
clear data_ICA;
% load options structure used to obtain the solution
data_options = load([input_dir, 'options.mat']);
options = data_options.options;
clear data_options;
if options.scen.unit_output == "mm"
    % convert S in m, so that the slip potency is in m^3 (SI units)
    ICA.S = ICA.S * mm2m;
end
% load the indices of the ICs selected for inversion
data_ind_comps = load([input_dir, 'ind_comps.mat']);
ind_comps_calculated = data_ind_comps.ind_comps;
clear data_ind_comps;
% load results of misfit from leave-one-out or k-fold cross-validation 
data_misfit_comps = load([input_dir, 'misfit_comps.mat']);
misfit_comps_calculated = data_misfit_comps.misfit_comps;
clear data_misfit_comps

% load plot parameters
run(files.plot_parameters);

% plot Power Spectral Density and ICs selection
options.inversion.select_comps.frequency_analysis.f_psd_st = 0.1;
options.inversion.select_comps.frequency_analysis.cs_psd_st = 0.9;
options.inversion.select_comps.common_mode_stddist = 3;
[fig_PSD, fig_cumPSD, ind_comps] = ...
    plot_PSD_selection(ICA, 'ICA', options, dirs);
close all;

[intersected_comps, ind_in_ind_comps_calculated] = ...
    intersect(ind_comps_calculated,ind_comps);

misfit_comps = misfit_comps_calculated(:,ind_in_ind_comps_calculated);

% plot ICs (space and time) for all ICs
fig_handle = plot_comp(Xd,ICA,[],[],options,[],[]);
close all;

% plot cross-validation curves (misfit vs smoothing parameter)
fig_misfit_comps = plot_select_smoothing(misfit_comps, options);
close all;

% load fault geometry and Green's functions
[fault, G] = load_geom(Xd, options, dirs);

% select best smoothing parameter for inversion
[min_misfit_comps,ind_sigma0_comps] = min(misfit_comps);

% invert selected ICs with selected smoothing parameter
[m, Cm] = invert_comps(ICA, ind_comps, fault, G, ind_sigma0_comps, ...
    options);

% slip potency components
n_ICs2invert = length(ind_comps);
n_patches = length(fault.area);
% slip potency components
mp = repmat(fault.area * km2m^2,2,1) .* m;
% co-variance matrix for slip potency components
Cmp = Cm;
parfor i=1:n_ICs2invert
    Cmp{i} = repmat(fault.area * km2m^2,2,1) .* Cm{i};
end

% create slip potency model
slip_potency = create_slip_model(mp, Cmp, ICA, fault, options, ind_comps);



options.plot.offsets = [];

options.plot.hypocenters.lon   = [];
options.plot.hypocenters.lat   = [];
options.plot.hypocenters.depth = [];
options.plot.hypocenters.mag   = [];
  
options.plot.borders = {'Canada', 'United States'};

options.plot.plot_err = 1;

options.plot.save_fig = 1;
options.plot.dir_save = [dirs.dir_raid,dirs.dir_scen,dirs.dir_case,...
    '/figures/'];

options.plot.FontSize = 16;

options.plot.comp.scale_factor = 0.02;
options.plot.comp.names        = 0;

%% PSD
options.plot.PSD.FontSize       = 32;
options.plot.PSD.FontSize_small = 18;
options.plot.PSD.LineWidth      = 2;
options.plot.PSD.time_unit      = 'yr';
options.plot.PSD.n_subplots_x   = 4;
options.plot.PSD.n_subplots_y   = 2;
options.plot.PSD.fig_pos        = [10 10];
options.plot.PSD.fig_size       = [1600 1200];
options.plot.PSD.output_name    = 'PSD';

options.plot.select_comps.fig_pos     = [10 10];
options.plot.select_comps.fig_size    = [600 600];
options.plot.select_comps.output_name = 'cumulativePSD';
options.plot.select_comps.output_name = 'Fig_select_comps';

options.plot.min_slip_rate = 0;

options.plot.map_slip_rate.FontSize     = 16;
options.plot.map_slip_rate.n_lats_edges = 101;
options.plot.map_slip_rate.fig_pos      = [10 10];
options.plot.map_slip_rate.fig_size     = [1800 800];
options.plot.map_slip_rate.colormap     = 'b2r'; % 'b2r' or 'hot'
options.plot.map_slip_rate.output_name  = 'Fig_slipratemap';


options.plot.video.FontSize         = 16;
options.plot.video.tremors_size     = 10;
options.plot.video.xlims            = [-128.5, -121.0];
options.plot.video.fig_pos          = [10 10];
options.plot.video.fig_size         = [418 800];
options.plot.video.field2plot_style = 'scatter'; % 'scatter' or 'fill'
options.plot.video.scatter_size     = 10;
options.plot.video.colormap         = 'b2r'; % 'b2r' or 'hot'
options.plot.video.colormap_minmax  = true;
options.plot.video.quality          = 100;
options.plot.video.frame_rate       = 10;
options.plot.video.quiver           = false;
options.plot.video.frac_quiver      = 0.02;
options.plot.video.output_name      = 'slip_rate';

options.plot.dynamics.FontSize         = 16;
options.plot.dynamics.fig_pos          = [10 10];
options.plot.dynamics.fig_size         = [800 800];
options.plot.dynamics.scatter_size_min = 10;
options.plot.dynamics.scatter_size_max = 50;
options.plot.dynamics.n_points_plot_d1theta = 7;
options.plot.dynamics.colormap_color   = 'hot';
options.plot.dynamics.colormap_minmax  = true;
options.plot.dynamics.fillvar          = 'max'; % 'max', 'min', 'sum'


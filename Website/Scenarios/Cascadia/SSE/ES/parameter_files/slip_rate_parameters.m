
options.slip_rate.rake = 90;

% sample rate (1/dt) in decimal years (dt is one day, i.e. dt = 1/365.25)
options.slip_rate.fs = 365.25;

% passband frequency in units of 1/yr (i.e., all periods lower than fpass
% yr ~ fpass*365.25 days are filtered)
options.slip_rate.fpass = 5;

% options.slip_rate.method can be:
% 'tvdiff', 'slopemovwin', 'slopemovwincen'
options.slip_rate.method = 'slopemovwincen';
switch options.slip_rate.method
    case 'tvdiff'
        options.slip_rate.tvdiff.iter  = 40;
        options.slip_rate.tvdiff.alph  = 0.1;
        options.slip_rate.tvdiff.scale = 'small';
        options.slip_rate.tvdiff.ep    = 1e-9;
        options.slip_rate.tvdiff.dx    = 'dt';
        options.slip_rate.tvdiff.plotflag = 0;
        options.slip_rate.tvdiff.diagflag = 0;
    case 'slopemovwin'
        options.slip_rate.slopemovwin.win_length = 28;
    case 'slopemovwincen'
        options.slip_rate.slopemovwincen.win_length = 7;
    otherwise

end
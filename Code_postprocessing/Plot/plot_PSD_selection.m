function [fig_PSD, fig_cumPSD, ind_comps] = plot_PSD_selection(...
    decomp, decomp_name, options, dirs)

FontSize        = options.plot.PSD.FontSize;
FontSize_small  = options.plot.PSD.FontSize_small;
LineWidth       = options.plot.PSD.LineWidth;
time_unit       = options.plot.PSD.time_unit;
n_subplots_x    = options.plot.PSD.n_subplots_x;
n_subplots_y    = options.plot.PSD.n_subplots_y;
fig_pos_PSD     = options.plot.PSD.fig_pos;
fig_size_PSD    = options.plot.PSD.fig_size;
output_name_PSD = options.plot.PSD.output_name;

n_freq        = options.inversion.select_comps.frequency_analysis.n_freq;
n_sample      = options.inversion.select_comps.frequency_analysis.n_sample;
flag_log_plot = ...
    options.inversion.select_comps.frequency_analysis.flag_log_plot;
output_name_select_comps = options.plot.select_comps.output_name;
fig_pos_select_comps     = options.plot.select_comps.fig_pos;
fig_size_select_comps    = options.plot.select_comps.fig_size;

f_skip = options.inversion.select_comps.frequency_analysis.f_skip;
sigmaf_skip=options.inversion.select_comps.frequency_analysis.sigmaf_skip;
f_threshold=options.inversion.select_comps.frequency_analysis.f_threshold;
cs_psd1 = options.inversion.select_comps.frequency_analysis.cs_psd1;
cs_psd2 = options.inversion.select_comps.frequency_analysis.cs_psd2;

f_psd_st = options.inversion.select_comps.frequency_analysis.f_psd_st;
cs_psd_st =options.inversion.select_comps.frequency_analysis.cs_psd_st;


n_comp = size(decomp.S,1);

if ~exist('time_unit','var')
    time_unit = 'yr';
end
if ~exist('decomp_name','var')
    decomp_name = 'decomp';
end
if ~exist('n_freq','var')
    n_freq = 50;
end
if ~exist('n_sample','var')
    n_sample = 2;
end
if ~exist('flag_log_plot','var')
    flag_log_plot = 0;
end

[f,PSDY,fit_log_PSDY] = frequency_analysis(decomp.V, options, ...
    decomp.timeline(1), decomp.timeline(end), time_unit, dirs, ...
    decomp_name, n_freq, n_sample, flag_log_plot);

n_subplots = n_subplots_x*n_subplots_y;
n_per_subplot = floor(options.scen.N / n_subplots);
n_extra = mod(options.scen.N, n_subplots);

fig_PSD = figure('Position', [fig_pos_PSD(1) fig_pos_PSD(2) ...
    fig_size_PSD(1) fig_size_PSD(2)]);

t = tiledlayout(n_subplots_y,n_subplots_x,'TileSpacing','Compact');
title(t,'ICs'' PSD','FontSize',FontSize)
xlabel(t,['Frequency (1/',time_unit,')'],'FontSize',FontSize)
ylabel(t,'PSD','FontSize',FontSize)
i = 0;
k = 0;
for iy=1:n_subplots_y
    for ix=1:n_subplots_x
        i = i+1;
        nexttile
        if i<=n_extra
            n_per_subplot_i = n_per_subplot+1;
        else
            n_per_subplot_i = n_per_subplot;
        end
        %subplot(n_subplots_y,n_subplots_x,i)
        %subaxis(n_subplots_y, n_subplots_x, i, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
        for j=1:n_per_subplot_i
            k = k + 1;
            if j<=2
                plot(f,PSDY{k}(1:numel(f)),'.-','LineWidth',LineWidth,...
                    'DisplayName', ['IC_{',num2str(k),'}']);
            else
                plot(f,PSDY{k}(1:numel(f)),'.--','LineWidth',LineWidth,...
                    'DisplayName', ['IC_{',num2str(k),'}']);
            end
            if j==1
                hold on;
                grid on;
            end
        end
        legend('show', 'Location', 'northeast')
        a = get(gca,'YTickLabel');
        set(gca,'YTickLabel',a,'FontSize',FontSize_small)
    end
end

if options.plot.save_fig == 1
    print([options.plot.dir_save,output_name_PSD],'-dpng','-r300');
end

%%%
% n_freq_f = numel(f);
% psd = cell(n_comp,1);
% for i=1:n_comp
%     psd{i} = PSDY{i}(1:n_freq_f);
% end
% 
% cs = zeros(n_comp,n_freq_f);
% scs = zeros(n_comp,1);
% for i=1:n_comp
%     psd_norm = psd{i} / sum(psd{i});
%     cs(i,:) = cumsum(psd_norm);
%     scs(i) = sum(cs(i,:));
% end
% 
% [f_thr,ind_f_thr] = min(abs(f-f_threshold));
% ind_comps = find(cs(:,ind_f_thr)>=cs_psd1 & cs(:,ind_f_thr)<=cs_psd2);
%%%

%%%
% f_skip = [1,2]; %options.inversion.select_comps.frequency_analysis.f_skip; % = [1,2];
% sigmaf_skip = [0.05, 0.05]; %options.inversion.select_comps.frequency_analysis.sigmaf_skip; %= [0.05, 0.05];
% f_threshold = 1.95; %options.inversion.select_comps.frequency_analysis.f_threshold;
% cs_psd1 = 0.5; %options.inversion.select_comps.frequency_analysis.cs_psd1;
% cs_psd2 = 0.95; %options.inversion.select_comps.frequency_analysis.cs_psd2;

n_freq_f = numel(f);
psd = cell(n_comp,1);
for i=1:n_comp
    psd{i} = PSDY{i}(1:n_freq_f);
end


f_low = f_skip-sigmaf_skip;
f_high = f_skip+sigmaf_skip;
ind_f2keep = true(size(f));
for i=1:size(f_low,2)
    ind_f2rm_tmp = f>=f_low(i) & f<=f_high(i);
    ind_f2keep(ind_f2rm_tmp) = false;
end
n_freq_f2keep = sum(ind_f2keep);
cs = zeros(n_comp,n_freq_f2keep);
for i=1:n_comp
    psd_norm = psd{i}(ind_f2keep) / sum(psd{i});
    cs(i,:) = cumsum(psd_norm);
end

[f_thr,ind_f_thr] = min(abs(f(ind_f2keep)-f_threshold));
ind_comps_cs_lt = find(cs(:,ind_f_thr)>=cs_psd1 & cs(:,ind_f_thr)<=cs_psd2);

[f_thr_st,ind_f_thr_st] = min(abs(f(ind_f2keep)-f_psd_st));
ind_comps_cs_st = find(((cs(:,ind_f_thr_st) - cs(:,1)) / ...
                        (f(ind_f_thr_st) - f(1))) < cs_psd_st );

ind_comps_cs = intersect(ind_comps_cs_st, ind_comps_cs_lt);

ind_comps2rm   = [];
for i=1:length(ind_comps_cs)
    j = ind_comps_cs(i);
    [psd_max, ind_f_max] = max(psd{j});
    f_max = f(ind_f_max);
    for k=1:length(f_skip)
        if (f_max >= f_skip(k)-sigmaf_skip(k)) && ...
                (f_max <= f_skip(k)+sigmaf_skip(k))
            ind_comps2rm = [ind_comps2rm, i];
        end
    end
end

ind_comps_complex = ind_comps_cs;
ind_comps_complex(ind_comps2rm) = [];


[ind_comps_common] = select_common_modes(decomp, options);
ind_comps_common_removed = intersect(ind_comps_complex, ind_comps_common);


if ~isempty(ind_comps_common_removed)
    ind_comps = ...
            ind_comps_complex(ind_comps_common_removed~=ind_comps_complex);
else
    ind_comps = ind_comps_complex;
end

%%%

fig_cumPSD = figure('Position', ...
    [fig_pos_select_comps(1) fig_pos_select_comps(2) ...
    fig_size_select_comps(1) fig_size_select_comps(2)]);
ax = gca();
hold on;
grid on;
plot(f(1:ind_f_thr), cs(:,1:ind_f_thr), '-k');
plot(f(1:ind_f_thr), cs(ind_comps,1:ind_f_thr), '-r', ...
        'LineWidth', LineWidth);
plot(f(1:ind_f_thr), cs(ind_comps_cs(ind_comps2rm),1:ind_f_thr), '--b', ...
        'LineWidth', LineWidth);
xlabel(['Frequency (1/', time_unit,')']);
ylabel('Cumulative PSD');
ax.FontSize = FontSize_small;
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'FontSize',FontSize_small)

if options.plot.save_fig == 1
    print([options.plot.dir_save,output_name_select_comps],'-dpng','-r300');
end

end
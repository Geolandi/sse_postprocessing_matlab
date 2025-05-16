function fig = plot_slip_rate_map(timeline, slip_rate, fault, tremors, ...
    options)

n_samples = size(slip_rate,2);

slip_rate_threshold_plot = options.plot.min_slip_rate;

n_lats_edges = options.plot.map_slip_rate.n_lats_edges;
colormap_color = options.plot.map_slip_rate.colormap;
FontSize = options.plot.map_slip_rate.FontSize;
fig_pos = options.plot.map_slip_rate.fig_pos;
fig_size = options.plot.map_slip_rate.fig_size;
output_name = options.plot.map_slip_rate.output_name;

n_lats = n_lats_edges - 1;
lat_min = min(min(fault.lat));
lat_max = max(max(fault.lat));
lats_edges = linspace(lat_min, lat_max, n_lats_edges);
lats = lats_edges(1:end-1) + (lats_edges(2:end)-lats_edges(1:end-1))/2;
slip_rate_lats_min = zeros(n_lats,n_samples);
slip_rate_lats_max = zeros(n_lats,n_samples);
%slip_rate_lats_sum = zeros(n_lats,n_samples);
for i=1:n_lats
    ind_lats = find(fault.llh(:,2)>lats_edges(i) & ...
        fault.llh(:,2)<=lats_edges(i+1));
    slip_rate_lats_min(i,:) = min(slip_rate(ind_lats,:));
    slip_rate_lats_max(i,:) = max(slip_rate(ind_lats,:));
    %slip_rate_lats_sum(i,:) = sum(slip_rate(ind_lats,:));
end

fig = figure('Position', [fig_pos(1) fig_pos(2) fig_size(1) fig_size(2)]);
hold on;
[T,L] = meshgrid(timeline, lats);

slip_rate_lats_minmax = slip_rate_lats_max;
switch colormap_color
    case 'hot'
        inds_rm = find(slip_rate_lats_max<slip_rate_threshold_plot);
        slip_rate_lats_minmax(inds_rm) = 0;
    case 'b2r'
        inds = find(slip_rate_lats_max<abs(slip_rate_lats_min));
        slip_rate_lats_minmax(inds) = slip_rate_lats_min(inds);
    otherwise

end
surf(T, L, zeros(n_lats, n_samples), slip_rate_lats_minmax,...
    'EdgeColor', 'None', 'FaceColor', 'flat')
if isempty(tremors)
    title('SSEs', 'FontSize', FontSize);
else
    scatter(tremors.timeline, tremors.lat, 1, 'k', 'filled')
    title('SSEs and tremors', 'FontSize', FontSize);
end
xlabel('Time (yr)', 'FontSize', FontSize);
ylabel('Lat (°)', 'FontSize', FontSize);
switch colormap_color
    case 'hot'
        colormap(flip(hot));
    case 'b2r'
        % colormap(b2r(-max(max(slip_rate_lats_minmax)), ...
        %     max(max(slip_rate_lats_minmax))))
        colormap(b2r(min(min(slip_rate_lats_minmax)), ...
            max(max(slip_rate_lats_minmax))))
    otherwise

end
cb = colorbar();
ylabel(cb,'Slip rate max & min (mm/yr)','FontSize',FontSize);
xlim([timeline(1), timeline(end)]);
ylim([lats_edges(1), lats_edges(end)]);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontSize',FontSize);
if options.plot.save_fig == 1
    print([options.plot.dir_save,output_name],'-dpng','-r300');
end

end
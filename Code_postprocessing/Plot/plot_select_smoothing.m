function fig = plot_select_smoothing(misfit_comps, options)

FontSize = options.plot.FontSize;
fig = figure;
title('Smoothing selection')
for i=1:size(misfit_comps,2)
    loglog(options.inversion.sigma, misfit_comps(:,i),'.-k');
    if i==1
        hold on;
        xlabel('\sigma_{m_0}', 'FontSize', FontSize);
        ylabel('Misfit spatial distribution (non-dimensional)', ...
            'FontSize', FontSize);
    end
end
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontSize',FontSize);

if options.plot.save_fig == 1
    print([options.plot.dir_save,'Fig_select_smoothing'],'-dpng','-r300');
end

end
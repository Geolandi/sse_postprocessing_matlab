function [f,PSDY,fit_log_PSDY] = frequency_analysis(V,options,first_epoch,last_epoch,time_unit,dirs,decomp_name,n_freq,n_sample,flag_log_plot)

% This function plot the Power Density Spectrum of the temporal
% eigenvectors of the PCA or the sources of the ICA.
% Input
% V: matrix containing as columns the temporal eigenvectors for a PCA or
% the sources for a ICA
% first_epoch: first epoch in the timeline
% last_epoch: last epoch in the timeline
% time_unit: time unit of the timeline
% n_freq: parameter that determines the number of frequency points to use
% for the computation of the Power Density Spectrum
% n_sample: additional parameter that allows to use a more refined sampling
% in the frequency domain

dir_scen = dirs.dir_scen;
dir_case = dirs.dir_case;

[n_time,n_comp] = size(V);
if n_freq < 2
    n_freq = 2;
    disp('Setting n_freq to 2: the other values will be symmetric.');
end
if n_freq > n_time
    n_freq = round(n_time/2);
    disp('Setting n_freq to round(n_time/2): you need at least two points.');
end
T = last_epoch - first_epoch;
switch time_unit
    case 'yr'
        f = 1/T*(0:round(365.25*T/n_freq));
        f = f/n_sample;
    case 'day'
        f = 1/T*(0:round(T/n_freq));
        f = f/n_sample;
    otherwise
        error('The time_unit is not supported');
end

Y = cell(1,n_comp);
PSDY = cell(1,n_comp);
fit_log_PSDY = cell(1,n_comp);
for ii=1:n_comp
    Y{ii} = fft(V(:,ii),n_time*n_sample);
    PSDY{ii} = Y{ii}.*conj(Y{ii})/numel(Y{ii});
    figure; hold on;
    grid on;
    title(['Power Spectral Density - ',decomp_name]);
    if flag_log_plot == 0
        plot(f,PSDY{ii}(1:numel(f)),'.-')
        xlabel(['Frequency (1/',time_unit,')']);
        ylabel(['Power Spectral Density for V_{',num2str(ii),'}']);
        if options.plot.save_fig == 1
            %saveas(gcf,[dir_scen,'figures/PSD/PSD_freq_',decomp_name,'_V',num2str(ii)],'fig');
            %saveas(gcf,[dir_scen,'figures/PSD/PSD_freq_',decomp_name,'_V',num2str(ii)],'png');
            %saveas(gcf,[dir_scen,dir_case,'figures/PSD/PSD_freq_',decomp_name,'_V',num2str(ii)],'pdf');
            saveas(gcf,[options.plot.dir_save,'PSD/PSD_freq_',decomp_name,'_V',num2str(ii)],'fig');
            %print([options.plot.dir_save,'PSD/PSD_freq_',decomp_name,'_V',num2str(ii)],'-dpng','-r300');
            %export_fig('-painters','-r600','-q101',[dir_scen,'figures/PSD/PSD_freq_',decomp_name,'_V',num2str(ii)],'pdf');
        end
    else
        plot(log(f),log(PSDY{ii}(1:numel(f))),'.-')
        xlabel(['Log frequency (with frequency in 1/',time_unit,')']);
        ylabel(['Log Power Spectral Density for V_{',num2str(ii),'}']);
        ind_notinf = find(isinf(log(f))==0);
        fit_log_PSDY{ii} = fit(log(f(ind_notinf))',log(PSDY{ii}(1:numel(f(ind_notinf)))),'poly1');
        plot(log(f(ind_notinf)),fit_log_PSDY{ii}.p1*log(f(ind_notinf)) + fit_log_PSDY{ii}.p2,'.-r');
        %fit_log_PSDY{ii} = fit(log(f(3:end))',log(PSDY{ii}(3:numel(f))),'poly1');
        %plot(log(f(3:end)),fit_log_PSDY{ii}.p1*log(f(3:end)) + fit_log_PSDY{ii}.p2,'.-r');
        axis equal
        if options.plot.save_fig == 1
            %saveas(gcf,[dir_scen,'figures/PSD/LogPSD_Logfreq_',decomp_name,'_V',num2str(ii)],'fig');
            %saveas(gcf,[dir_scen,'figures/PSD/LogPSD_Logfreq_',decomp_name,'_V',num2str(ii)],'png');
            %saveas(gcf,[dir_scen,dir_case,'figures/PSD/LogPSD_Logfreq_',decomp_name,'_V',num2str(ii)],'pdf');
            saveas(gcf,[options.plot.dir_save,'PSD/LogPSD_Logfreq_',decomp_name,'_V',num2str(ii)],'fig');
            %print([options.plot.dir_save,'PSD/LogPSD_Logfreq_',decomp_name,'_V',num2str(ii)],'-dpng','-r300');
            %export_fig('-painters','-r600','-q101',[dir_scen,'figures/PSD/LogPSD_Logfreq_',decomp_name,'_V',num2str(ii)],'pdf');
        end
    end
    
    
    
end

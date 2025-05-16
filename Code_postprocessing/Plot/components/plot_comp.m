function fig_handle = ...
    plot_comp(X,decomp,model,fig_handle,options,faults,seismicity)

% plot_comp.m creates a number of figures equal to the number of components
% specified in the structure decomp. It plots for each component the
% temporal evolution of the source, and the correpsonding map associated
% with the mixing matrix. The source is scaled so that is oscillates in the
% range [0,1]. After such a rescaling, the spatial distribution is also
% rescaled, incorporating also the weight that gives the unit of
% measurement.
% -------------------------------------------------------------------------
% INPUT
% decomp: Structure that must contain the following fields:
%           U: Matrix MxN. Spatial distribution
%           S: Matrix NxN. Weights
%           V: Matrix TxN. Temporal functions
%           var_U: Matrix MxN. Variance on the spatial distribution
%           var_V: Matrix TxN. Variance on the temporal functions
%           type: String specifying the type of data used for the
%                 decomposition.
%           llh: Matrix Mx3. Longitude, latitude, and height for the
%                spatial distributions.
%           timeline: Vector 1xT. Epochs for the temporal functions.
% fig_handle: Figure handle variable containing the figure handles for all
%             already existing figures.
% options: Optional. Structure that must contain the following fields:
%          plot: Structure that must contain the following fields:
%                offsets: Vector containing the epochs at which a green
%                         line will be drawn in the temporal function plot.
%                borders: Cell array containing the borders to plot )empty
%                         for no borders)
% faults: Optional. Structure that must contain the following field:
%         faults_model: Cell array with N elements, one for every
%                       components in the decomposition. If a component is
%                       not flagged as to be inverted, no geometry is
%                       loaded for it and the faults_model will be empty.
%                       Otherwise it contains a cell of size Nfaults. This
%                       cell contains Nfaults matrices, each of them with
%                       a specific fault geometry in the PCAIM format.
%         traces: Structure loaded by the load_faults_traces_kml.m file.
% seismicity: Optional. Structure containing the seismicity. It will be
%             plotted on map as circles. For the temporal part, the
%             normalized cumulative number of events is plotted.
% -------------------------------------------------------------------------
% OUTPUT
% fig_handle: Updated figure handle variable.
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 6 Nov 2017
% Jet Propulsion Laboratory, California Institute of Technology

scale_factor = options.plot.comp.scale_factor;

% Set font size
fontsize = 14;

if nargin<3
    flag_model = 0;
    model = [];
else
    if ~isempty(model)
        flag_model = 1;
    else
        flag_model = 0;
    end
end
if nargin<4
    fig_handle = [];
end
if nargin<5
    options = [];
    borders_list = [];
else
    borders_list = options.plot.borders;
end
if nargin<6
    faults.faults_model = [];
    faults.traces       = [];
end
if nargin<7
    seismicity = [];
end

% if ~isempty(options.faults.kmlfile2plot)
%     faults.traces = load_faults_traces_kml([],options);
% else
%     faults.traces = [];
% end

Nborders = numel(borders_list);
% Extract usefuel variables from the structure decomp
[U,S,V,var_U,var_V,types,llh,timeline,N,M,T] = ...
    extract_variables_from_decomp(decomp,options);

if options.scen.unit_output == "mm"
    mm2m = 1e-3;
    S = S / mm2m;
else
    S = S;
end

if isempty(faults)
    faults.faults_model = cell(1,N);
    faults.traces       = [];
    faults.whichfaults  = cell(1,N);
else
    faults.whichfaults  = options.faults.whichfaults;
end

%type = types{1}(1:4);
type = types{1};
if strcmp(type(1:4),'GPS3')
    type(5:end) = [];
end

signV = sign(V(end,:) - V(1,:));
V = V.*repmat(signV,T,1);
U = U.*repmat(signV,M,1);

% Define the epoch of the offsets to plot
toffsets = options.plot.offsets;
% Define the number of offsets
Noffsets = numel(toffsets);

if options.plot.comp.names==1
    % Define the names of the processed time series in X
    names_direction = X.name;
    M = numel(names_direction);
    names_repeated = cell(M,1);
    for jj=1:M
        try
            names_repeated{jj} = names_direction{jj}(1:4);
        catch
            names_repeated{jj} = num2str(names_direction{jj});
        end
    end
    names = unique(names_repeated);
end
% Determine the max and min of the temporal functions
switch decomp.decmode
    case 'T-mode'
        maxtfunc_notnorm = max(V);
        mintfunc_notnorm = min(V);
        
        tfunc_notnorm = V;
        sfunc_notnorm = U*S;
        var_tfunc_notnorm = var_V;
        var_sfunc_notnorm = var_U*(S.^2);
        if flag_model==1
            sfunc_modeled_notnorm = model.data_modeled*S;
        end
    case 'S-mode'
        A = U*S;
        maxtfunc_notnorm = max(A);
        mintfunc_notnorm = min(A);
        
        tfunc_notnorm = A;
        sfunc_notnorm = V;
        var_tfunc_notnorm = var_U*(S.^2);
        var_sfunc_notnorm = var_V;
        if flag_model==1
            sfunc_modeled_notnorm = model.data_modeled*S;
        end
        clear A
    otherwise
        
end

% Find x_lim and y_lim for the map plot
x_lim_map = [min(llh(1:3:end,1)) - ...
    (max(llh(1:3:end,1)) - min(llh(1:3:end,1)))/20, ...
    max(llh(1:3:end,1)) + ...
    (max(llh(1:3:end,1)) - min(llh(1:3:end,1)))/20];
y_lim_map = [min(llh(1:3:end,2)) - ...
    (max(llh(1:3:end,2)) - min(llh(1:3:end,2)))/20,...
    max(llh(1:3:end,2)) + ...
    (max(llh(1:3:end,2)) - min(llh(1:3:end,2)))/20];

% For every component...
for ii=1:N
    % Determine the scaling factor for the temporal function
    facttfunc = maxtfunc_notnorm(ii)-mintfunc_notnorm(ii);
    % Set the temporal function such that it is between 0 and 1
    tfunc(:,ii) = (tfunc_notnorm(:,ii)-mintfunc_notnorm(ii))./facttfunc;
    var_tfunc(:,ii) = var_tfunc_notnorm(:,ii)./(facttfunc^2);
    err_tfunc(:,ii) = sqrt(var_tfunc(:,ii));
    % Rescale the spatial function to maintain the product with the
    % temporal function
    sfunc(:,ii) = sfunc_notnorm(:,ii)*facttfunc;
    var_sfunc(:,ii) = var_sfunc_notnorm(:,ii).*(facttfunc^2);
    err_sfunc(:,ii) = sqrt(var_sfunc(:,ii));
    if flag_model==1
        sfunc_modeled(:,ii) = sfunc_modeled_notnorm(:,ii)*facttfunc;
    end
    % Create a figure
    figure;
    fig_handle = [fig_handle; gcf];
    hold on; grid on;
    % Subdivide it in four, and fill the first row with the temporal source
    ax1 = subplot(4,1,1); hold on; grid on;
    if options.plot.plot_err==1
        errorbar(ax1,timeline,tfunc(:,ii),err_tfunc(:,ii),'.r');
    else
        plot(ax1,timeline,tfunc(:,ii),'.r');
    end
    ylim([-0.5,1.5]);
    y_lim = get(gca,'YLim');
    for kk=1:Noffsets
        line(ax1,[toffsets(kk), toffsets(kk)],y_lim,'Color','b','LineWidth',4,'LineStyle','--');
    end
    %title(['Component number ',num2str(ii)]);
    xlabel('Time (yr)');
    set(ax1,'xaxisLocation','top')
    switch X.decmode
        case 'T-mode'
            ylabel(['source_{',num2str(ii),'}']);
        case 'S-mode'
            ylabel(['Mixing matrix A_{',num2str(ii),'}']);
        otherwise
            
    end
    
    % Fill the remaining three quarters of the figure with the
    % spatial map
    ax2 = subplot(4,1,[2,4]); hold on; grid on;
    xlim(x_lim_map);
    ylim(y_lim_map);
    for kk=1:Nborders
        % Plot the specified borders
        if ~isempty(borders_list{kk})
            borders(borders_list{kk},'b','nomap');
        end
    end
    xlabel('Longitude ({\circ})','Interpreter','tex');
    ylabel('Latitude ({\circ})','Interpreter','tex');
    
    if ~isempty(faults.whichfaults{ii})
        flags_plot_fault.flag_llh   = 1;
        flags_plot_fault.flag_white = 1;
        plot_faults(faults,model,ax2,flags_plot_fault,options);
        view(2)
    end
    
    ax3 = axes;
    xlim(x_lim_map);
    ylim(y_lim_map);
    if ~isempty(faults.whichfaults{ii})
        if ~isempty(faults.traces)
            plot_fault_traces(faults,options,ax3);
            view(2)
        end
    end
    
    ax4 = axes;
    xlim(x_lim_map);
    ylim(y_lim_map);
    if ~isempty(seismicity)
        if ~isempty(faults.whichfaults{ii})
            scatter3(ax4,seismicity.lon,seismicity.lat,...
                seismicity.depth,10 + ...
                50*log10((min(seismicity.mag) + ...
                seismicity.mag)/max(seismicity.mag)),...
                seismicity.timeline,'filled')
            view(2)
        end
    end
    
    ax5 = axes;
    xlim(x_lim_map);
    ylim(y_lim_map);
    if options.plot.comp.names==1
        llh_stn = llh(1:3:end,:);
        for jj=1:M/3
            text(ax5,llh_stn(jj,1),llh_stn(jj,2),names{jj});
        end
    end
    
    switch type
        case 'GPS3'
            
            ax8 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            % EXTRA VERTICAL LEGEND
            scatter(ax8,-114.8,33.4,500,0,'Filled','MarkerEdgeColor','k');
            scatter(ax8,-114.8,33.3,200,0,'Filled','MarkerEdgeColor','k');
            text(ax8,-115,33.4,'Data','FontSize',fontsize);
            text(ax8,-115.05,33.3,'Model','FontSize',fontsize);
            
            hscatter = scatter(ax8,llh(1:3:end,1),llh(1:3:end,2),500,...
                sfunc(3:3:end,ii),'Filled','MarkerEdgeColor','k');
            if flag_model==1
                if ~isempty(faults.whichfaults{ii})
                    hold on;
                    scatter(ax8,llh(1:3:end,1),llh(1:3:end,2),200,...
                        sfunc_modeled(3:3:end,ii),'Filled',...
                        'MarkerEdgeColor','k');
                end
            end
            view(2)
            
            ax6 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            llh_legend = [min(llh(1:3:end,1)) + ...
                (max(llh(1:3:end,1)) - min(llh(1:3:end,1)))/20,...
                min(llh(1:3:end,2)) + ...
                (max(llh(1:3:end,2)) - min(llh(1:3:end,2)))/20];
            
            vector_field_legend = ...
                10*ceil(0.1*0.5*max(sqrt(sfunc(1:3:end,ii).^2 + sfunc(2:3:end,ii).^2)));
            vector_field = scale_factor*[sfunc(1:3:end,ii),sfunc(2:3:end,ii); ...
                vector_field_legend, 0];
            llh_vector_field = [llh(1:3:end,1:2); llh_legend];
            quiver(ax6,llh_vector_field(:,1),llh_vector_field(:,2),...
                vector_field(:,1),vector_field(:,2),'r','LineWidth',4,...
                'AutoScale','on');
            if flag_model==1
                hold on;
                llh_legend_modeled = [llh_legend(1),llh_legend(2)-0.1];
                llh_vector_field_modeled = [llh(1:3:end,1:2); ...
                    llh_legend_modeled];
                vector_field_modeled = scale_factor*...
                    [sfunc_modeled(1:3:end,ii),sfunc_modeled(2:3:end,ii); ...
                    vector_field_legend, 0];
                quiver(ax6,llh_vector_field_modeled(:,1),...
                    llh_vector_field_modeled(:,2),...
                    vector_field_modeled(:,1),...
                    vector_field_modeled(:,2),'b',...
                    'LineWidth',4,'AutoScale','on');
            end
            
            
            ax7 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            text(ax7,llh_legend(1),llh_legend(2)+...
                2*(max(llh(1:3:end,1)) - min(llh(1:3:end,1)))/20,...
                'Horizontal','FontSize',fontsize);
            text(ax7,llh_legend(1),llh_legend(2)+...
                (max(llh(1:3:end,1)) - min(llh(1:3:end,1)))/20,...
                [num2str(vector_field_legend),' mm'],'FontSize',fontsize);
            text(ax7,llh_legend(1)+...
                2.5*(max(llh(1:3:end,2)) - min(llh(1:3:end,2)))/50,...
                llh_legend(2),'Data','FontSize',fontsize);
            if flag_model==1
                text(ax7,llh_legend_modeled(1)+...
                    2.5*(max(llh(1:3:end,2)) - min(llh(1:3:end,2)))/50,...
                    llh_legend_modeled(2),'Model','FontSize',fontsize);
            end
            
            
            ax9 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            if ~isempty(faults.whichfaults{ii})
                if ~isempty(options.plot.hypocenters.lat)
                    scatter3(options.plot.hypocenters.lon,...
                        options.plot.hypocenters.lat,...
                        -options.plot.hypocenters.depth,...
                        2000 + 150*log10((min(options.plot.hypocenters.mag) + ...
                        options.plot.hypocenters.mag)/...
                        max(options.plot.hypocenters.mag)),'p','filled','g',...
                        'MarkerEdgeColor','k');
                end
            end
            view(2)
            
            linkaxes([ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9])
            % Hide the top axes
            ax3.Visible = 'off';
            ax3.XTick = [];
            ax3.YTick = [];
            ax4.Visible = 'off';
            ax4.XTick = [];
            ax4.YTick = [];
            ax5.Visible = 'off';
            ax5.XTick = [];
            ax5.YTick = [];
            ax6.Visible = 'off';
            ax6.XTick = [];
            ax6.YTick = [];
            ax7.Visible = 'off';
            ax7.XTick = [];
            ax7.YTick = [];
            ax8.Visible = 'off';
            ax8.XTick = [];
            ax8.YTick = [];
            ax9.Visible = 'off';
            ax9.XTick = [];
            ax9.YTick = [];
            
            % Give each one its own colormap
            if flag_model==1
                colormap(ax2,flipud(hot));
                colormap(ax8,b2r(...
                    min([min(sfunc(3:3:end,ii)),min(sfunc_modeled(3:3:end,ii))]),...
                    max([max(sfunc(3:3:end,ii)),max(sfunc_modeled(3:3:end,ii))])));
            else
                %colormap(ax2,'jet')
                colormap(ax8,b2r(min(sfunc(3:3:end,ii)),max(sfunc(3:3:end,ii))));
            end
            colormap(ax4,'cool')
            
            % Then add colorbars and get everything lined up
            set(ax1,'Position',[0.1 0.725 .8 .2],'FontSize',fontsize);
            set([ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],'Position',[0.1 0.1 .8 .6],'FontSize',fontsize);
            
            if ~isempty(faults.whichfaults{ii})
                cb2 = colorbar(ax2,'Position',[.75 .12 .03 .25]);
                if ~isempty(seismicity)
                    cb4 = colorbar(ax4,'Position',[.75 .45 .03 .25]);
                end
            end
            cb8 = colorbar(ax8,'Position',[.2 .3 .03 .25]);
            
            hscatter.addprop('X');
            hscatter.X = X;
            hscatter.addprop('decomp');
            hscatter.decomp = decomp;
            hscatter.addprop('options');
            hscatter.options = options;
            hscatter.addprop('nn');
            hscatter.nn = ii;
        
        case 'c-well'
            
            ax8 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            hscatter = scatter(ax8,llh(:,1),llh(:,2),500,...
                sfunc(:,ii),'Filled','MarkerEdgeColor','k');
            if flag_model==1
                if ~isempty(faults.whichfaults{ii})
                    hold on;
                    scatter(ax8,llh(:,1),llh(:,2),200,...
                        sfunc_modeled(:,ii),'Filled',...
                        'MarkerEdgeColor','k');
                end
            end
            view(2)
            
            ax6 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            
            ax7 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            
            ax9 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            if ~isempty(faults.whichfaults{ii})
                if ~isempty(options.plot.hypocenters.lat)
                    scatter3(options.plot.hypocenters.lon,...
                        options.plot.hypocenters.lat,...
                        -options.plot.hypocenters.depth,...
                        2000 + 150*log10((min(options.plot.hypocenters.mag) + ...
                        options.plot.hypocenters.mag)/...
                        max(options.plot.hypocenters.mag)),'p','filled','g',...
                        'MarkerEdgeColor','k');
                end
            end
            view(2)
            
            linkaxes([ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9])
            % Hide the top axes
            ax3.Visible = 'off';
            ax3.XTick = [];
            ax3.YTick = [];
            ax4.Visible = 'off';
            ax4.XTick = [];
            ax4.YTick = [];
            ax5.Visible = 'off';
            ax5.XTick = [];
            ax5.YTick = [];
            ax6.Visible = 'off';
            ax6.XTick = [];
            ax6.YTick = [];
            ax7.Visible = 'off';
            ax7.XTick = [];
            ax7.YTick = [];
            ax8.Visible = 'off';
            ax8.XTick = [];
            ax8.YTick = [];
            ax9.Visible = 'off';
            ax9.XTick = [];
            ax9.YTick = [];
            
            % Give each one its own colormap
            if flag_model==1
                colormap(ax2,flipud(hot));
                colormap(ax8,b2r(...
                    min([min(sfunc(:,ii)),min(sfunc_modeled(:,ii))]),...
                    max([max(sfunc(:,ii)),max(sfunc_modeled(:,ii))])));
            else
                %colormap(ax2,'jet')
                colormap(ax8,b2r(min(sfunc(:,ii)),max(sfunc(:,ii))));
            end
            colormap(ax4,'cool')
            
            % Then add colorbars and get everything lined up
            set(ax1,'Position',[0.1 0.725 .8 .2],'FontSize',fontsize);
            %set([ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],'Position',[0.1 0.1 .6 .6],'FontSize',fontsize);
            set([ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],'Position',[0.1 0.1 .8 .6],'FontSize',fontsize);
            
            if ~isempty(faults.whichfaults{ii})
                cb2 = colorbar(ax2,'Position',[.75 .12 .03 .25]);
                if ~isempty(seismicity)
                    cb4 = colorbar(ax4,'Position',[.75 .45 .03 .25]);
                end
            end
            %cb8 = colorbar(ax8,'Position',[.85 .25 .01 .25]);
            cb8 = colorbar(ax8,'Position',[.75 .43 .03 .25]);
            
            hscatter.addprop('X');
            hscatter.X = X;
            hscatter.addprop('decomp');
            hscatter.decomp = decomp;
            hscatter.addprop('options');
            hscatter.options = options;
            hscatter.addprop('nn');
            hscatter.nn = ii;
        case 'InSARLOS'
            ax8 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            
            ax9 = axes;
            xlim(x_lim_map);
            ylim(y_lim_map);
            if ~isempty(faults.whichfaults{ii})
                if ~isempty(options.plot.hypocenters.lat)
                    scatter3(options.plot.hypocenters.lon,...
                        options.plot.hypocenters.lat,...
                        -options.plot.hypocenters.depth,...
                        2000 + 150*log10((min(options.plot.hypocenters.mag) + ...
                        options.plot.hypocenters.mag)/...
                        max(options.plot.hypocenters.mag)),'p','filled','g',...
                        'MarkerEdgeColor','k');
                end
            end
            view(2)
            
            linkaxes([ax2,ax3,ax4,ax5,ax8,ax9])
            % Hide the top axes
            ax3.Visible = 'off';
            ax3.XTick = [];
            ax3.YTick = [];
            ax4.Visible = 'off';
            ax4.XTick = [];
            ax4.YTick = [];
            ax5.Visible = 'off';
            ax5.XTick = [];
            ax5.YTick = [];
            ax8.Visible = 'off';
            ax8.XTick = [];
            ax8.YTick = [];
            ax9.Visible = 'off';
            ax9.XTick = [];
            ax9.YTick = [];
            
            if options.flags.flag_downsample == 0
                hscatter = scatter(ax8,llh(:,1),llh(:,2),40,...
                    sfunc(:,ii),'s','Filled');
            else
                lon_patch = zeros(size(llh,1),4);
                lat_patch = zeros(size(llh,1),4);
                for jj=1:size(llh,1)
                    lon_min = min(X.polys{jj}(:,2));
                    lon_max = max(X.polys{jj}(:,2));
                    lat_min = min(X.polys{jj}(:,3));
                    lat_max = max(X.polys{jj}(:,3));
                    lon_patch(jj,:) = [lon_min, lon_max, lon_max, lon_min];
                    lat_patch(jj,:) = [lat_min, lat_min, lat_max, lat_max];
                end
                hscatter = patch(ax8,lon_patch',lat_patch',sfunc(:,ii),'EdgeColor','None');
            end
            drawnow;
            if flag_model==1
                if ~isempty(faults.whichfaults{ii})
                    hold on;
                    scatter(ax8,llh(:,1),llh(:,2),20,...
                        sfunc_modeled(:,ii),'s','Filled');
                end
            end
            view(2)
            if max(abs(sfunc(:,ii)))==0
                colormap(ax8,b2r(-eps,eps));
                caxis(ax8,[-eps,eps]);
            else
                colormap(ax8,b2r(-max(abs(sfunc(:,ii))),max(abs(sfunc(:,ii)))));
                caxis(ax8,[-max(abs(sfunc(:,ii))),max(abs(sfunc(:,ii)))]);
            end
            cb8 = colorbar(ax8,'Position',[.91 .43 .03 .25]);
            hscatter.addprop('X');
            hscatter.X = X;
            hscatter.addprop('decomp');
            hscatter.decomp = decomp;
            hscatter.addprop('options');
            hscatter.options = options;
            hscatter.addprop('nn');
            hscatter.nn = ii;
            
            % Then add colorbars and get everything lined up
            set(ax1,'Position',[0.1 0.725 .8 .2],'FontSize',fontsize);
            set([ax2,ax3,ax4,ax5,ax8,ax9],'Position',[0.1 0.1 .8 .6],'FontSize',fontsize);
        
        otherwise
        
    end
    if ~isempty(faults.whichfaults{ii})
        if flag_model==1
            set(get(cb2,'ylabel'),'string','Slip (mm)','FontSize',fontsize);
        else
            set(get(cb2,'ylabel'),'string','Depth (km)');
        end
        if ~isempty(seismicity)
            set(get(cb4,'ylabel'),'string','Time (yr)');
        end
    end
    switch type
        case 'cGPS3'
            set(get(cb8,'ylabel'),'string','Vertical (mm)','FontSize',fontsize);
        case 'c-well'
            set(get(cb8,'ylabel'),'string','Water level (mm)','FontSize',fontsize);
        case 'InSARLOS'
            set(get(cb8,'ylabel'),'string','LOS (mm)','FontSize',fontsize);
        otherwise
            
    end
    
    %daspect([1 1 1]);
    set(hscatter,'buttonDownFcn',@plot_ts_stn_click);
    set(fig_handle(end),'Resize','off')
    screen_size = get(groot,'Screensize');
    l_fig = min(0.8*[screen_size(3),screen_size(4)]);
    set(fig_handle(end), 'Position', [100, 100, l_fig, l_fig]);
    pause(2);
    drawnow;
    if options.plot.save_fig==1
        %saveas(gcf,[options.plot.dir_save,'components/comp_',num2str(ii),'.fig'],'fig');
        print([options.plot.dir_save,'components/comp_',num2str(ii)],'-dpng','-r300');
    end
end

end
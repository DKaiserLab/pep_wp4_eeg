function plot_corr_sig(cfg, d)

% evaluate input
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 20; end
if ~isfield(cfg, 'p_type'); cfg.p_type = 'p_vals'; end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.05, 0.08]; end
if ~isfield(cfg, 'xlim'); cfg.xlim = [-200, 500]; end
cfg.plot_preds = 1;


% create figure and layout
fig = figure;
set(gcf, 'Units','centimeters', 'Position',[2 2 80 20]);
tiledlayout(1, length(cfg.frequencies))

% graphic parameters
set(0, 'DefaultTextFontSize', cfg.FontSize);
set(0, 'DefaultAxesFontSize', cfg.FontSize);
set(0, 'DefaultTextFontName', cfg.FontName);
set(0, 'DefaultAxesFontName', cfg.FontName);
set(gcf, 'color', [1 1 1]);


% loop over frequency bands
for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};


    % prep time points and variables
    timepoints = d.pairRep.(frqBand).included_time;
    timeseries = d.pairRep.(frqBand).all_time;
    x = timeseries(timepoints)*1000;

    nexttile
    title(['IS-RSA - ', frqBand], 'FontSize', cfg.FontSize+10, 'FontWeight', 'bold');
    hold on;

    % prep time points and variables
    plot_counter = 1;
    g = gobjects(1, length(cfg.RDM_to_partial_out));

    for var = 1:cfg.plot_preds
        % smooth data
        % plot mean

        y = d.r_vals.ISC_RSA.(frqBand).r_obs_mean;
        if strcmp(frqBand, 'full')
            y = smoothdata(y, 'movmean', cfg.smoothing_window);
        end

        % compute signifikant time points
        sigMat = d.stats.ISC_RSA.(frqBand).(cfg.p_type);
        if isempty(sigMat)
            sigMat = false(length(x), 1);
        else
            sigMat = sigMat < 0.05;
        end 

        % get color
        if startsWith(cfg.RDM_to_partial_out{var}, 'typical')
            clr = [1, 0, 1];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'control')
            clr = [.7, .7, .7];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')
            clr = [.4, .9, 1];
        end

        
        % plot emprical distribution [5%, 95%] percentile
        e = d.stats.ISC_RSA.(frqBand).ci;
        boundedline(x, zeros(length(x)), abs(e), 'cmap', [1, 0, 1]);
        
        
        % plot mean
        g(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 5, ...
            'DisplayName', regexprep(cfg.RDM_to_partial_out{var}, "_.*", ""));

        % mark signifikant time points
        pos = max(y) + 0.01;
        plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
            'color', clr, 'marker' ,'O', 'MarkerFaceColor', clr ,'MarkerSize', 7, 'LineStyle','none');

        plot_counter = plot_counter + 1;
    end

    if cfg.partial_cor
        ylabel('Partial correlation');
    else
        ylabel([cfg.correlation_type, ' correlation [r]', newline]);
    end

    ylim(cfg.ylim);
    xlim(cfg.xlim);
    yline(0, '-', 'LineWidth', 2);
    xline(0, '--', 'LineWidth', 2);
    xlabel('Time (ms)');
    set(gca, 'box', 'off');
    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);
    hold off;

end
end
function plot_corr_sig(cfg, d, pval)

if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6; end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.1, 0.1]; end

fig = figure;
tiledlayout(1,2);   % 1 Zeile, 2 Spalten

timepoints = d.pairRep.included_time;
timeseries = d.pairRep.all_time;
x = timeseries(timepoints);
x_lim = ([min(x)-min(x)-0.01, max(x)+0.01]);
ntp_before_stim = length(d.pairRep.included_time) - length(d.stats.testedTime);

for category = cfg.categories
    nexttile;
    hold on;
    category = char(category);
    legend_labels = cell(1,3);
    plot_counter = 1;
    p = gobjects(1, length(cfg.RDM_to_partial_out)*length(cfg.categories));
    for var = 1:length(cfg.RDM_to_partial_out)

        % smooth data
        y = d.resMat.partial_cor.(category)(var, :);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

        % compute significant time points
        if isfield(pval, category)
            sigmat = pval.(category)(var,:) <= 0.05;
        else
            sigmat = pval.FDR.(category)(var,:);
        end

        sigMat = [false(1,ntp_before_stim), sigmat];

        % get color
        if startsWith(cfg.RDM_to_partial_out{var}, 'typical')
            clr = [1, 0, 1];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'control')
            clr = [.7, .7, .7];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')
            clr = [.4, .9, 1];
        end

        % plot line
        p(plot_counter) = plot(x, y, 'color', clr, 'LineWidth', 3);

        % mark signifikant time points
        pos = 0.1 - var*0.01;
        plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
            'color', clr, 'marker' ,'*', 'MarkerSize', 5, 'LineStyle','none');

        legend_labels{plot_counter} = [strrep(cfg.RDM_to_partial_out{var}, '_', '-'), ' ', category];
        plot_counter = plot_counter+1;
    end


    xlim(x_lim);
    ylim(cfg.ylim);
    yline(0, '--')
    xline(0, '--')
    xlabel('time');
    legend([p(1:plot_counter-1)], legend_labels, 'Location', 'best');
    set(gca, 'box', 'off');

end % category

sgtitle('Time-resolved correlation ISC-RDM with predictors', 'FontSize', 18);
hold off;
save_plot(fig, 'corr_sig_pvalFDR_early', cfg.figPath)

% 
fig = figure;
hold on;
legend_labels = cell(1,3);
plot_counter = 1;
p = gobjects(1, length(cfg.RDM_to_partial_out));

for var = 1:numel(cfg.RDM_to_partial_out)
    % smooth data
    % plot mean
    y1 = d.resMat.partial_cor.bathroom(var, :);
    y2 = d.resMat.partial_cor.kitchen(var, :);
    y = mean([y1; y2]);
    y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

    % compute signifikant time points
    if isfield(pval, category)
        sigmat = pval.all(var,:) <= 0.05;
    else
        sigmat = pval.FDR.all(var,:);
    end

    sigMat = [false(1, ntp_before_stim), sigmat];

    % get color
    if startsWith(cfg.RDM_to_partial_out{var}, 'typical')
        clr = [1, 0, 1];
    elseif startsWith(cfg.RDM_to_partial_out{var}, 'control')
        clr = [.7, .7, .7];
    elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')
        clr = [.4, .9, 1];
    end

    % plot mean
    p(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 3);

    % mark signifikant time points
    pos = 0.1 - var*0.01;
    plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
        'color', clr, 'marker' ,'*', 'MarkerSize', 5, 'LineStyle','none');

    legend_labels{plot_counter} = [strrep(cfg.RDM_to_partial_out{var}, '_', '-'),' mean'];
    plot_counter = plot_counter + 1;
end

if cfg.partial_cor
    ylabel(['Partial correlation [r]', newline]);
else
    ylabel([cfg.correlation_type, ' correlation [r]', newline]);
end

set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold');

ylim(cfg.ylim);
xlim(x_lim);
yline(0, '--');
xline(0, '--');
xlabel('time');
legend([p(1:plot_counter-1)], legend_labels, 'Location', 'best');
set(gca, 'box', 'off');

sgtitle('Time-resolved correlation ISC-RDM with predictors', 'FontSize', 18);
hold off;

save_plot(fig, 'mean_corr_sig_pvalFDR_early', cfg.figPath)

end
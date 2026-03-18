function p_val_fdr = plot_decoding_results(cfg, d)


if ~isfield(cfg, 'chance_level'); cfg.chance_level = 0.5; end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'xlim'); cfg.xlim = [-200, 500]; end

% Initialize data storage
p_val = struct();
t_stat = struct();
is_significant = struct();
p_val_fdr = struct();

% initialize figure
figure;
set(gcf, 'Units','centimeters', 'Position',[2 2 80 20]);
if length(cfg.frequencies) > 1
    tl = tiledlayout(1, length(cfg.frequencies));
    title(tl,'Mean pairwise decoding accuracy', 'FontWeight', 'bold', 'FontSize', cfg.FontSize +1)
end

% loop over frequency bands
for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};

    % get timepoints
    x = d.pairRep.(frqBand).all_time(d.pairRep.(frqBand).included_time)*1000;
    n_timepoints = length(d.pairRep.(frqBand).included_time);

    for i = 1:length(cfg.categories)
        cat = cfg.categories{i};

        % Initialize data storage
        p_vals_cat = zeros(1, n_timepoints);
        t_vals_cat = zeros(1, n_timepoints);

        % compute t-test for each time point
        for t = 1:n_timepoints
            [~, p, ~, stats] = ttest(d.meanAcc.(cat).(frqBand)(t,:), cfg.chance_level, 'tail', 'right');
            p_vals_cat(t) = p;
            t_vals_cat(t) = stats.tstat;
        end

        % FDR
        [is_significant.(cat).(frqBand), ~, ~, adj_p] = fdr_bh(p_vals_cat);

        % save
        p_val.(cat).(frqBand) = p_vals_cat;
        t_stat.(cat).(frqBand) = t_vals_cat;
        p_val_fdr.(cat).(frqBand) = adj_p;

    end

    % create next subplot
    if length(cfg.frequencies) > 1
        nexttile
    end
    set(gcf, 'Units','centimeters', 'Position',[2 2 40 25]);
    hold on;

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    set(gcf, 'color', [1 1 1]);

    % preallocate
    h = gobjects(1, length(cfg.categories));

    colors = lines(length(cfg.categories));
    pos = (max([mean(d.meanAcc.kitchen.(frqBand),2);...
        mean(d.meanAcc.bathroom.(frqBand), 2)]) + 0.016) * 100; % position to mark significant time points

    for i=1:length(cfg.categories)

        category = cfg.categories{i};
        sem = (std(d.meanAcc.(category).(frqBand), 0, 2)' ./ sqrt(cfg.n))* 100;
        mean_acc = mean(d.meanAcc.(category).(frqBand), 2)' * 100;
        upper = mean_acc + sem;
        lower = mean_acc - sem;

        fill([x, fliplr(x)], ...
            [upper, fliplr(lower)], colors(i, :), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h(i) = plot(x, mean_acc, 'Color', colors(i, :), 'LineWidth', 5, 'DisplayName', category);

        % mark significant time points
        for t = 1:n_timepoints
            if is_significant.(category).(frqBand)(t) == 1
                plot(x(t), pos-0.2*i,...
                    'marker' ,'O', 'MarkerFaceColor', colors(i, :) , 'Color', colors(i, :), 'MarkerSize', 7)
            end
        end
    end

    xlim(cfg.xlim)
    ylim([49, pos+0.6])
    xline(0, '--', 'LineWidth', 2);
    yline(50, '--', 'LineWidth', 2);

    % add labels and legend
    xlabel('Time (ms)');
    ylabel('Decoding accuracy (%)');
    lgd = legend(h, 'Location', 'none');
    set(lgd, 'Position', [0.7 0.7 0.3 0.1]);
    legend('boxoff');
    if length(cfg.frequencies) > 1
        title(frqBand);
    else
        title('Mean pairwise decoding accuracy')
    end

    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);
    set(gca, 'box', 'off');

    % save figure
    %     if cfg.saving
%         save_plot(fig, 'pairwise-decoding-with-sig-tp', cfg.figPath)
%     end

end
end

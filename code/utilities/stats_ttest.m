function stats_ttest(cfg, d)

if ~isfield(cfg, 'chance_level'); cfg.chance_level = 0.5; end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end

% Initialize data storage
n_timepoints = size(d.pairRep.included_time, 2);
p_val = struct();
t_stat = struct();
is_significant = struct();
p_val_fdr = struct();

for i = 1:length(cfg.categories)
    cat = cfg.categories{i};
    
    % Initialize data storage
    p_vals_cat = zeros(1, n_timepoints);
    t_vals_cat = zeros(1, n_timepoints);

    % compute t-test for each time point
    for t = 1:n_timepoints
        [~, p, ~, stats] = ttest(d.meanAcc.(cat)(t,:), cfg.chance_level, 'tail', 'right');
        p_vals_cat(t) = p;
        t_vals_cat(t) = stats.tstat;
    end

    % FDR
    [is_significant.(cat), ~, ~, adj_p] = fdr_bh(p_vals_cat, 0.050);

    % save
    p_val.(cat) = p_vals_cat;
    t_stat.(cat) = t_vals_cat;
    p_val_fdr.(cat) = adj_p;

end

% stats.p_val = p_val;
% stats.t_stat = t_stat;
% stats.p_val_fdr = p_val_fdr;

if cfg.plotting

    fig = figure;
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
    x = d.pairRep.all_time(d.pairRep.included_time);
    y_min = min([mean(d.meanAcc.kitchen,2); mean(d.meanAcc.bathroom,2)]) - 0.002; % position to mark significant time points

    for i=1:length(cfg.categories)
        
        cat = cfg.categories{i};
        sem = std(d.meanAcc.(cat), 0, 2)' ./ sqrt(cfg.n);
        mean_acc = mean(d.meanAcc.(cat),2)';
        upper = mean_acc + sem;
        lower = mean_acc - sem;

        fill([x, fliplr(x)], ... 
            [upper, fliplr(lower)], colors(i, :), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h(i) = plot(x, mean_acc, 'Color', colors(i, :), 'LineWidth', 2);

        % mark significant time points
        plot(x(is_significant.(cat)), repmat(y_min-0.002*strcmp(cat, 'bathroom'),1,sum(is_significant.(cat))), '*', 'Color', colors(i, :), 'MarkerSize', 3);
    end

    xlim([min(x)-min(x)-0.01, max(x)+0.01])
    xline(0, '--');
    yline(0.5, '--');

    % add labels and legend
    xlabel('Time (s)');
    ylabel('Mean decoding accuracy');
    legend(h, cfg.categories, 'Location', 'best');
    title(sprintf('Mean pairwise decoding accuracy\n with FDR-corrected significant time points'));

    % save figure
    if cfg.saving
        save_plot(fig, 'pairwise-decoding-with-sig-tp', cfg.figPath)
    end
end

end

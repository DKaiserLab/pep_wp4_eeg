function stats_ttest_stimuli(cfg, d)
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05; end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6; end

% Initialize data storage
n_timepoints = length(d.pairRep.included_time);
p_val = struct();
t_stat = struct();
is_significant = struct();
p_val_fdr = struct();

for l=1:length(cfg.RDM_to_partial_out)
    layer = cfg.RDM_to_partial_out{l};
    for i = 1:length(cfg.categories)
        cat = cfg.categories{i};

        % Initialize data storage
        p_vals_cat = struct;
        t_vals_cat = struct;
        p_vals_cat.(layer) = zeros(1, n_timepoints);
        t_vals_cat.(layer) = zeros(1, n_timepoints);

        % compute t-test for each time point
        for iTp = 1:n_timepoints
            [~, p, ~, stats] = ttest(d.stimuli.all.(layer).(cat)(:,iTp), 0, 'tail', 'right');
            p_vals_cat.(layer)(iTp) = p;
            t_vals_cat.(layer)(iTp) = stats.tstat;
        end

        % FDR
        [is_significant.(layer).(cat), ~, ~, adj_p] = fdr_bh(p_vals_cat.(layer), cfg.alpha);

        % save
        p_val.(layer) = struct;
        t_stat.(layer) = struct;
        p_val.(layer).(cat) = p_vals_cat;
        t_stat.(layer).(cat) = t_vals_cat;
        p_val_fdr.(layer) = struct;
        p_val_fdr.(layer).(cat) = adj_p;

    end % cat
end % layer

% mean
for l=1:length(cfg.RDM_to_partial_out)
    layer = cfg.RDM_to_partial_out{l};
    y_all = (d.stimuli.all.(layer).bathroom + d.stimuli.all.(layer).kitchen)/2;
    n_timepoints = size(y_all, 2);
    p_mean = zeros(1, n_timepoints);
    t_mean = zeros(1, n_timepoints);

    for iTp = 1:n_timepoints
        [~, p_mean(iTp), ~, stats] = ttest(y_all(:, iTp), 0, 'tail', 'right');
        t_mean(iTp) = stats.tstat;
    end

    % fdr
    is_significant_mean.(layer) = fdr_bh(p_mean, cfg.alpha);
end


%% plotting
if cfg.plotting
    % evtl. noch sem hinzufügen

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    set(gcf, 'color', [1 1 1]);

    fig = figure;
    tiledlayout(length(cfg.RDM_to_partial_out),1);

    x=d.pairRep.all_time(d.pairRep.included_time);
    colors = lines(length(cfg.categories));
    
    for i=1:length(cfg.RDM_to_partial_out)
        nexttile; hold on;
        
        %preallocate
        lower = nan(length(cfg.categories),n_timepoints);
        upper = nan(length(cfg.categories),n_timepoints);
        p = gobjects(1,length(cfg.categories));

        layer = cfg.RDM_to_partial_out{i};
        for cat=1:length(cfg.categories)
            category = cfg.categories{cat};

            % compute and plot mean
            mean_cat = mean(d.stimuli.all.(layer).(category),1);
            p(cat)=plot(x, mean_cat, 'LineWidth', 1, 'Color', colors(cat,:));

            % plot standard error of the mean
            sem = std(d.stimuli.all.(layer).(category),0, 1)/sqrt(cfg.n);
            upper(cat,:) = mean_cat + sem;
            lower(cat,:) = mean_cat - sem;
            fill([x, fliplr(x)], ...
                [upper(cat,:), fliplr(lower(cat,:))], colors(cat,:),...
                'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
        end
       
        % mark significant time points
        for cat=1:length(cfg.categories)
            category = cfg.categories{cat};
            h_idx = find(is_significant.(layer).(category));
            plot(x(h_idx), repmat(max(max(upper))+0.005*cat, 1, length(h_idx)), '*', 'Color', colors(cat,:), 'MarkerSize', 3, 'HandleVisibility', 'off');
        end

        title(strrep(layer, '_', '-'));
        ylim([-0.05, 0.1]);
        xlim([min(x)-0.01, max(x)+0.01]);
        yline(0,'--');
        xline(0,'--');
        xlabel('time (s)');
        ylabel(cfg.correlation_type)
        set(gca, 'box', 'off');

    end
    sgtitle('Correlation of time point RDMs with stimuli RDM');
    lg = legend(p, cfg.categories);
    lg.Layout.Tile = 'east';
    % save plot
    if cfg.saving
        save_plot(fig, cfg.save_name, cfg.figPath);
    end

    %%
    figure;
    hold on;
    p = gobjects(1, length(cfg.RDM_to_partial_out));  % Länge an die Schleife anpassen
    legend_labels = cell(1, length(cfg.RDM_to_partial_out));
    for i = 1:length(cfg.RDM_to_partial_out)

        layer = cfg.RDM_to_partial_out{i};
        % smooth data
        % plot mean
        y1 = mean(d.stimuli.all.(layer).bathroom,1);
        y2 = mean(d.stimuli.all.(layer).kitchen,1);
        y = mean([y1; y2]);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

        % get color
        if endsWith(cfg.RDM_to_partial_out{i}, 'early') %strcmp(cfg.RDM_to_partial_out{var}, 'typical_late')
            clr = [1, 0, 1];
        elseif endsWith(cfg.RDM_to_partial_out{i}, 'mid') %strcmp(cfg.RDM_to_partial_out{var}, 'control_late')
            clr = [.7, .7, .7];
        elseif endsWith(cfg.RDM_to_partial_out{i}, 'late')%strcmp(cfg.RDM_to_partial_out{var}, 'photos_late')
            clr = [.4, .9, 1];
        end
        % plot mean
        p(i) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 3);
        legend_labels{i} = [strrep(cfg.RDM_to_partial_out{i}, '_', '-'),' mean'];
        h_idx = find(is_significant_mean.(layer));
        plot(x(h_idx), repmat(0.06+0.01*i, 1, length(h_idx)), '*', 'Color', clr, 'MarkerSize', 8, 'HandleVisibility','off');
        
        ylabel(cfg.correlation_type);
        set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold');
        
    end

    legend(p, legend_labels, 'Location','northeastoutside');
    ylim([-0.05, 0.1]);
    xlim([min(x)-0.01, max(x)+0.01]);
    yline(0,'--', 'HandleVisibility','off');
    xline(0,'--', 'HandleVisibility','off');
    xlabel('time (s)');
    ylabel(cfg.correlation_type)

    title('Correlation of time point RDMs with stimuli RDM');
    set(gca, 'box', 'off');

    % saving
    if cfg.saving
        fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], 'compare_tp_RDMs_to_stimuli_RDMs');
        fig_name = 'Time-resolved-corr-reference-RDMs-with-stimuli';
        save_plot(fig, fig_name, fig_path)
    end
end

end

function d = compare_tp_RDMs_to_predictor_RDMs(d, cfg)

% evaluate input
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson';end
if ~isfield(cfg, 'add_legend'); cfg.add_legend = true;end
if ~isfield(cfg, 'show_single_cate'); cfg.show_single_cate = false;end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true;end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_tp_RDMs_to_predictor_RDMs';end
if ~isfield(cfg, 'xaxis_labels'); cfg.xaxis_labels = true;end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
%if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.3, 0.5];end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6;end
cfg.plot_rdm = false;
cfg.ISC_type = 'pairRep';

% prepare figure
if cfg.plotting
    legend_labels = {};
    plot_counter = 1;
    fig = figure;
    hold on
end

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % init results table
    if cfg.partial_cor
        d.resMat.partial_cor.(category) = nan(numel(cfg.RDM_to_partial_out), numel(d.(cfg.ISC_type).included_time));
        d.resMat.partial_cor.model = cfg.RDM_to_partial_out;
    else
        d.resMat.cor.(category) = nan(1, numel(d.(cfg.ISC_type).included_time));
        d.resMat.cor.model = cfg.analysis_names;
    end
    % get canditate/predictor RDMs
    RDMs = struct;
    RDMs(1).name = d.ISC.([category,'_RDM']).(cfg.ISC_type).name;
    RDMs(1).color = d.ISC.([category,'_RDM']).(cfg.ISC_type).color;
    RDMs(1).RDM = d.ISC.([category,'_RDM']).(cfg.ISC_type)(1).RDM; % (:, :, 1);
    labels = {RDMs.name};
    [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);

    % loop through time points
    for iTp = 1:numel(d.(cfg.ISC_type).included_time)

        % get time point RDM
        RDMs(1).RDM = d.ISC.([category,'_RDM']).(cfg.ISC_type)(iTp).RDM; %(:, :, iTp);

        % partial correlation
        if cfg.partial_cor
            [~, rMat, ~, cfg] = partial_cor_RDM(cfg, RDMs);
        else
            [~, rMat, ~] = cor_RDM(RDMs,cfg);
        end

        if cfg.partial_cor
            d.resMat.partial_cor.(category)(:, iTp) = rMat(2:end, 1);
        else
            d.resMat.cor.(category)(iTp) = rMat;
        end

    end

    % plotting
    if cfg.plotting && cfg.partial_cor
        x = d.(cfg.ISC_type).included_time;

        for var = 1:numel(cfg.RDM_to_partial_out)

            % smooth data
            y = d.resMat.partial_cor.(category)(var, :);
            y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

            % get color
            if startsWith(cfg.RDM_to_partial_out{var}, 'typical') %strcmp(cfg.RDM_to_partial_out{var}, 'typical_late')
                clr = [1, 0, 1];
            elseif startsWith(cfg.RDM_to_partial_out{var}, 'control') %strcmp(cfg.RDM_to_partial_out{var}, 'control_late')
                clr = [.7, .7, .7];
            elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')%strcmp(cfg.RDM_to_partial_out{var}, 'photos_late')
                clr = [.4, .9, 1];
            end

            % plot line
            if strcmp(category, cfg.categories{1}) % bathroom
                p(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', '--', 'LineWidth', 1);
                yline(0, '--')
            else % kitchen
                p(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', ':', 'LineWidth', 1);
                yline(0, '--')
            end
            legend_labels{plot_counter} = [strrep(cfg.RDM_to_partial_out{var}, '_', '-'), ' ', category];
            plot_counter = plot_counter+1;
        end
    end
end

if cfg.plotting && cfg.partial_cor

    for var = 1:numel(cfg.RDM_to_partial_out)

        % smooth data
        % plot mean
        y1 = d.resMat.partial_cor.bathroom(var, :);
        y2 = d.resMat.partial_cor.kitchen(var, :);
        y = mean([y1; y2]);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

        % get color
        if startsWith(cfg.RDM_to_partial_out{var}, 'typical') %strcmp(cfg.RDM_to_partial_out{var}, 'typical_late')
            clr = [1, 0, 1];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'control') %strcmp(cfg.RDM_to_partial_out{var}, 'control_late')
            clr = [.7, .7, .7];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')%strcmp(cfg.RDM_to_partial_out{var}, 'photos_late')
            clr = [.4, .9, 1];
        end
        % plot mean
        p(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 3);

        legend_labels{plot_counter} = [strrep(cfg.RDM_to_partial_out{var}, '_', '-'),' mean'];

        if cfg.partial_cor
            ylabel(['Partial correlation [r]', newline]);
        else
            ylabel([cfg.correlation_type, ' correlation [r]', newline]);
        end

        set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold');
       
        plot_counter = plot_counter + 1;
    end

    % add legend
    if cfg.add_legend
        legend(p, legend_labels, 'Location','northeastoutside');
    end

    title('Time-resolved correlation reference RDMs with predictors')

    if isfield(cfg, 'ylim')
        ylim(cfg.ylim)
    end

    xlabel('time');
    set(gca, 'XTickLabel', round(d.pairRep.all_time(x([1, 10:10:60])), 2));
    set(gca, 'box', 'off');

    % saving
    if cfg.saving
        fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], 'compare_tp_RDMs_to_predictor_RDMs');
        fig_name = 'Time-resolved-corr-reference-RDMs-with-predictors';
        save_plot(fig, fig_name, fig_path)
    end
end
end
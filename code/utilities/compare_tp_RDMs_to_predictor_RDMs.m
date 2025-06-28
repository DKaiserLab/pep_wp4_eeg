function d = compare_tp_RDMs_to_predictor_RDMs(d, cfg)

% evaluate input
if ~isfield(cfg, 'RDM_to_partial_out'); cfg.RDM_to_partial_out = {'typical_late', 'control_late'}; end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson';end
if ~isfield(cfg, 'plot_rdm'); cfg.plot_rdm = false;end
if ~isfield(cfg, 'add_legend'); cfg.add_legend = true;end
if ~isfield(cfg, 'show_single_cate'); cfg.show_single_cate = false;end
if ~isfield(cfg, 'order_predictors'); cfg.order_predictors = false;end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = true;end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_roi_RDMs_to_predictor_RDMs';end
if ~isfield(cfg, 'xaxis_labels'); cfg.xaxis_labels = true;end
cfg.plot_rdm = false;
if ~isfield(cfg, 'plot_type'); cfg.plot_type = 'bar';end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.3, 0.5];end
if ~isfield(cfg, 'ISC_type'); cfg.ISC_type = 'pairRep';end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6;end

% prepare figure
if cfg.plotting
    figure;
    hold on
end

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % init results table
    d.resMat.(category) = nan(numel(cfg.RDM_to_partial_out), numel(d.(cfg.ISC_type).included_time));

    % get canditate/predictor RDMs
    RDMs = struct;
    RDMs(1).name = d.([category,'_RDM']).(cfg.ISC_type).name;
    RDMs(1).color = d.([category,'_RDM']).(cfg.ISC_type).color;
    RDMs(1).RDM = d.([category,'_RDM']).(cfg.ISC_type).RDM(:, :, 1);
    labels = {RDMs.name};
    [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);

    % loop through time points
    for iTp = 1:numel(d.(cfg.ISC_type).included_time)

        % get time point RDM
        RDMs(1).RDM = d.([category,'_RDM']).(cfg.ISC_type).RDM(:, :, iTp);

        % partial correlation
        if cfg.partial_cor
            [~, rMat, ~, cfg] = partial_cor_RDM(cfg, RDMs);
        else
            [~, rMat, ~] = cor_RDM(RDMs,cfg);
        end
        d.resMat.(category)(:, iTp) = rMat(2:end, 1);
    end

    % plotting
    if cfg.plotting
        x = d.(cfg.ISC_type).included_time;

        for var = 1:numel(cfg.RDM_to_partial_out)

            % smooth data
            y = d.resMat.(category)(var, :);
            y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

            % get color
            if strcmp(cfg.RDM_to_partial_out{var}, 'typical_late')
                clr = [1, 0, 1];
            elseif strcmp(cfg.RDM_to_partial_out{var}, 'control_late')
                clr = [.7, .7, .7];
            elseif strcmp(cfg.RDM_to_partial_out{var}, 'photos_late')
                clr = [.4, .9, 1];
            end

            % plot line
            if strcmp(category, cfg.categories{1})
                p(var) = plot(x, y, 'color', clr, 'LineStyle', '--', 'LineWidth', 1);
                yline(0, '--')
            else
                p(var) = plot(x, y, 'color', clr, 'LineStyle', ':', 'LineWidth', 1);
                yline(0, '--')
            end
        end
    end
end

if cfg.plotting

    for var = 1:numel(cfg.RDM_to_partial_out)

        % smooth data
        % plot mean
        y1 = d.resMat.bathroom(var, :);
        y2 = d.resMat.kitchen(var, :);
        y = mean([y1; y2]);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

        % get color
        if strcmp(cfg.RDM_to_partial_out{var}, 'typical_late')
            clr = [1, 0, 1];
        elseif strcmp(cfg.RDM_to_partial_out{var}, 'control_late')
            clr = [.7, .7, .7];
        elseif strcmp(cfg.RDM_to_partial_out{var}, 'photos_late')
            clr = [.4, .9, 1];
        end
        % plot mean
        p(end) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 3);

        if cfg.partial_cor
            ylabel(['Partial correlation [r]', newline]);
        else
            ylabel([cfg.correlation_type, ' correlation [r]', newline]);
        end

        title('Time-resolved correlation reference RDMs with predictors')
%         if isfield(cfg, 'ylim')
%             ylim(cfg.ylim)
%         end

        set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold')
        ax = gca;
        ax.Box = 'off';

        % add legend to last plot
        if cfg.add_legend
            legend(p, cfg.RDM_to_partial_out, 'Location','northeastoutside');
        end

        %     % saving
        %     fig_path = fullfile(pwd, 'figures', ['exp_', num2str(cfg.exp_num)], 'compare_roi_RDMs_to_predictor_RDMs');
        %     save_plot(cfg.save_name, fig_path)
    end
end
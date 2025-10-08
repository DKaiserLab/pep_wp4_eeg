function [d, stats] = getISCofRespresentation(cfg, d)

if ~isfield(cfg, 'plot_sem'); cfg.plot_sem = true;end
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'Spearman';end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_tp_RDMs_to_stimuli_RDMs';end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05;end
if ~isfield(cfg, 'n_permutations'); cfg.n_permutations = 5000;end
if ~isfield(cfg, 'permutationtest'); cfg.permutationtest = false;end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 1;end
if ~isfield(cfg, 'xlim'); cfg.xlim = [-200, 500]; end
cfg.plot_rdm = false;

ntimepoints= length(d.pairRep.included_time);
stats = struct;

% prepare permutation test
if cfg.permutationtest
    rng(1);
    signs_all = randi([0 1], cfg.n, cfg.n_permutations)*2 - 1;
    null_distrib = struct;
    stats.cluster = struct;

    for cate_num = 1:numel(cfg.categories)
        category = char(cfg.categories{cate_num});
        null_distrib.(category) = nan(ntimepoints, cfg.n_permutations);
        stats.cluster.(category) = struct;
    end

end

% loop through categories
for cate_num = 1:numel(cfg.categories)
    category = char(cfg.categories{cate_num});

    % loop through time points
    for tp=1:ntimepoints

        % preallocate RDM matrix
        RDMmat = nan(nchoosek(cfg.nTrials/2, 2), cfg.n);

        % loop through subjects
        for iSub = 1:length(cfg.subNums)
            subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
            subID2 = strrep(subID, '-', '');

            % make a matrix with vectorized RDMs
            rdm = squeeze(d.pairRep.(subID2).rdm(:,:,tp));
            rdm(eye(size(rdm)) == 1) = 0;

            % filter for category
            if strcmp(category, 'bathroom')
                rdm = rdm(1:cfg.nTrials/2, 1:cfg.nTrials/2);
            else
                rdm = rdm(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end);
            end

            RDMmat(:, iSub) = squareform(rdm);

        end

        % make IS-RDM
        [~, mat_out, ~] = make_RDM(RDMmat, cfg);
        if cfg.dissimilarity
            median_mat_out = 1 - mat_out;
        else
            median_mat_out = mat_out;
        end

        % take median and standard error
        median_mat_out(eye(size(median_mat_out)) == 1) = 0;
        medianISC = median(squareform(median_mat_out), 'omitnan'); % median VZ-permutieren, ganze zeile+spalte
        seISC = std(squareform(median_mat_out), 'omitnan') / sqrt(length(squareform(median_mat_out)));

        if cfg.permutationtest

            if isempty(gcp('nocreate'))
                parpool(8);
            end

            perm_median = nan(1, cfg.n_permutations);
            parfor perm = 1:cfg.n_permutations
                signs = signs_all(:,perm);
                signs_mat = diag(signs)
                perm_median_mat_out = signs_mat * median_mat_out * signs_mat;
                perm_medianISC = median(squareform(perm_median_mat_out), 'omitnan');
                perm_median(perm) = perm_medianISC;
            end
            null_distrib.(category)(tp,:) = perm_median;
        end
        
        % store in structure
        d.ISC.([category,'_RDM']).pairRep(tp).RDM = mat_out;
        d.ISC.([category,'_RDM']).pairRep(tp).color = [0, 0, 0];
        d.ISC.([category,'_RDM']).pairRep(tp).name = (num2str(tp));
        d.ISC.medianISC.(category).pairRep(tp) = medianISC;
        d.ISC.medianISC_SE.(category).pairRep(tp) = seISC;

    end % time point
    if cfg.permutationtest
        stats.cluster.(category) = stats_cluster(d.ISC.medianISC.(category).pairRep, null_distrib.(category), cfg, d);
    end
end % category




%% plotting
if cfg.plotting

    %% median ISC plot
    fig=figure;
    set(gcf, 'Units','centimeters', 'Position',[2 2 40 25]);
    hold on;
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    set(gcf, 'color', [1 1 1]);

    c = lines(numel(cfg.categories));
    x = d.pairRep.all_time(d.pairRep.included_time)*1000;
    h = gobjects(1, length(cfg.categories)); % +1
    max_ISC = round(max([d.ISC.medianISC.bathroom.pairRep, d.ISC.medianISC.kitchen.pairRep]), 2);
    min_ISC = round(min([d.ISC.medianISC.bathroom.pairRep, d.ISC.medianISC.kitchen.pairRep]), 2);
    
    for icat=1:length(cfg.categories)
        category = cfg.categories{icat};
        y = smoothdata(d.ISC.medianISC.(category).pairRep,'movmean', cfg.smoothing_window);
        se = d.ISC.medianISC_SE.(category).pairRep;

        % plot median ISC
        h(icat)=plot(x, y, 'color', c(icat, :), 'LineWidth', 5, 'DisplayName', category);
        if cfg.plot_sem
            x2 = [x, fliplr(x)];
            inBetween_bath = [y+ se, fliplr(y - se)];
            fill(x2, inBetween_bath, c(icat, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        
        if cfg.permutationtest
            sigMat = false(1, ntimepoints);
            if ~isempty(stats.cluster.(category).obs_cluster_stats)
                for p=1:length(stats.cluster.(category).pvals_cluster)
                    if stats.cluster.(category).pvals_cluster(p) <= cfg.alpha
                        idx = stats.cluster.(category).obs_clusters.PixelIdxList{p};
                        sigMat(idx)=1;
                    end
                end
            end

            pos = max_ISC + icat*0.002;
            plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
                'color', c(icat, :), 'marker' ,'O', 'MarkerFaceColor', c(icat, :) ,'MarkerSize', 7, 'LineStyle','none');
        end
        
    end

%     y1 = d.ISC.medianISC.bathroom.pairRep;
%     y2 = d.ISC.medianISC.kitchen.pairRep;
%     y = mean([y1; y2]);
%     y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);
%     h(icat+1)=plot(x, y, 'color', 'black', 'LineWidth', 2);

    set(gca, 'box', 'off');
    yline(0, '--', 'LineWidth', 2);
    xline(0, '--','LineWidth', 2);  % Mark stimulus onset
    xlabel('Time (ms)');
    ylabel('Spearman''s R');
    xlim(cfg.xlim)
    ylim([min_ISC, max_ISC+0.006]);
    
    %title({'Median inter-subject correlation', 'of neural representations'}, 'FontSize', cfg.FontSize + 10);
    legend(h);
    legend('boxoff');
    set(gca, 'box', 'off');
    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);
    if cfg.saving
        save_plot(fig, 'ISCofRepresentation', cfg.figPath);
    end

end

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

end
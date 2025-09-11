function cluster_stats = stats_cluster(r_obs, perm_corrs, cfg, d)

% evaluate input
if ~isfield(cfg,'alpha'); cfg.alpha = 0.05;end
if ~isfield(cfg,'min_cluster_size'); cfg.min_cluster_size = 3; end
if ~isfield(cfg,'plot_treshhold'); cfg.plot_treshhold = false; end
if size(r_obs,1) == 1
    r_obs = r_obs(:);
end

% define parameters
thresh_per_tp = prctile(perm_corrs, 100*(1-cfg.alpha), 2);
above_thresh = r_obs(:) > thresh_per_tp(:);
clusters = bwconncomp(above_thresh); % all time points above treshhold
valid_idx = cellfun(@numel, clusters.PixelIdxList) >= cfg.min_cluster_size; 
clusters.PixelIdxList = clusters.PixelIdxList(valid_idx); % just cluster with size >= min_cluster_size
clusters.NumObjects = sum(valid_idx); % number of significant clusters

obs_cluster_stats = zeros(1, clusters.NumObjects); % init stats
for c = 1:clusters.NumObjects
    obs_cluster_stats(c) = sum(r_obs(clusters.PixelIdxList{c})); % add correlations at given time points
end
perm_max_cluster = nan(1,cfg.n_permutations);

% find biggest cluster in permutations
for p = 1:cfg.n_permutations
    perm_stat = perm_corrs(:,p); % [timepoints x correlation], "correlation-distribution" over time
    
    % check if perm correlations > definded treshhold
    above = perm_stat > thresh_per_tp;
    clusts = bwconncomp(above);

    valid_idx = cellfun(@numel, clusts.PixelIdxList) >= cfg.min_cluster_size; 
    clusts.PixelIdxList = clusts.PixelIdxList(valid_idx); % just cluster with size >= min_cluster_size
    clusts.NumObjects = sum(valid_idx); % number of significant clusters

    % compute cluster-mass
    max_stat = 0;
    for c = 1:clusts.NumObjects
        s = sum(perm_stat(clusts.PixelIdxList{c}));
        if s > max_stat
            max_stat = s; % check if cluster-mass is higher then the previous one
        end
    end
    perm_max_cluster(p) = max_stat;
   
end

% compute p-values
pvals_cluster = nan(size(obs_cluster_stats));
for c = 1:length(obs_cluster_stats)
    pvals_cluster(c) = mean(perm_max_cluster >= obs_cluster_stats(c)); % check if perm-cluster-mass is higher then observed
end

% save output
cluster_stats = struct;
cluster_stats.obs_clusters = clusters;
cluster_stats.obs_cluster_stats = obs_cluster_stats;
cluster_stats.pvals_cluster = pvals_cluster;
cluster_stats.perm_max_cluster = perm_max_cluster;

if cfg.plot_treshhold
    figure; hold on;
    ntp_before_stim = length(d.pairRep.included_time) - length(d.stats.testedTime);
    x = d.pairRep.all_time(d.pairRep.included_time);
    x = x(ntp_before_stim+1:end);
    plot(x, r_obs, 'b-', 'LineWidth', 1.5);           % observed correlations
    plot(x, thresh_per_tp, 'r--', 'LineWidth', 1.2);  % treshhold for each timepoint

    % mark cluster
    colors = lines(clusters.NumObjects);
    for ci = 1:clusters.NumObjects
        idxs = clusters.PixelIdxList{ci} + ntp_before_stim+1;
        plot(idxs, r_obs(idxs), 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', colors(ci,:), ...
            'MarkerEdgeColor', 'k');
    end

    xlabel('Time points');
    ylabel('Observed r');
    title('cluster permutationtest');
    legend({'Observed r','Threshold'},'Location','best');
    grid on;
end

end

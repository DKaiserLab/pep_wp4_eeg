function tfce_ds=EEG_clusterstats_time_adapted(cfg,data)

    %requires:
    %cfg.time1
    %cfg.time2
    %cfg.nIter
    %cfg.h0
    
    %this adapted version assumes an input array (i.e. the data argument) 
    %that is subject num * time 1 * time 2

    ss=1:size(data,3);

    ds.sa.chunks=ss';
    ds.sa.targets=ones(size(ss))';

    allFeat=ones(length(ds_group),length(cfg.time2));

    ds.fa.train_time=allFeat.*[1:length(cfg.time1)]'; 
    ds.fa.train_time=ds.fa.train_time(:)'; 

    ds.fa.test_time=allFeat.*[1:length(cfg.time2)]; 
    ds.fa.test_time=ds.fa.test_time(:)';

    ds.a.fdim.values={1:length(cfg.time1),1:length(cfg.time2)};
    ds.a.fdim.labels={'train_time','test_time'};
    
    for s=ss
        ds.samples(s,:)=reshape(squeeze(data(:,:,s)),length(cfg.time1)*length(cfg.time2),1);
    end

    cluster_nbrhood=cosmo_cluster_neighborhood(ds);

    opt=struct();
    opt.niter=cfg.nIter;
    opt.h0_mean=cfg.h0;

    %this gives out a 1xcfg.time*cfg.time vector of z-values
    tfce_ds=cosmo_montecarlo_cluster_stat(ds,cluster_nbrhood,opt);
    
end
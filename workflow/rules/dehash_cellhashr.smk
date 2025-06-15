rule dehash_cellhashr:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample, "raw_feature_bc_matrix"),
        whitelist_path = lambda wc: os.path.join(RESULTS_DIR, "empty", wc.sample, "whitelist.txt"),
        htotosampletsv_path = lambda wc: os.path.join(HTOMAPPING_DIR, wc.sample, "hto_to_sample_mapping.tsv")
    output:
        barcode_metadata = os.path.join(RESULTS_DIR, "dehash", "{sample}", "barcode_metadata.csv"),
        metrics = os.path.join(RESULTS_DIR, "dehash", "{sample}", "metrics.csv"),
        output = os.path.join(RESULTS_DIR, "dehash", "{sample}", "output.csv"),
        session = os.path.join(RESULTS_DIR, "dehash", "{sample}", "session_info.txt"),
        whitelist = os.path.join(RESULTS_DIR, "dehash", "{sample}", "whitelist.txt")
    params:
        main = lambda wc: params_to_cli_args(config["dehash_cellhashr"]),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, "dehash", wc.sample),
    log:
        os.path.join(RESULTS_DIR, "dehash", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --counts_dir {input.counts_dir} \
            --htotosampletsv_path {input.htotosampletsv_path} \
            --out_dir {params.out_dir} \
            --whitelist_path {input.whitelist_path} \
            {params.main} &> {log}
        """

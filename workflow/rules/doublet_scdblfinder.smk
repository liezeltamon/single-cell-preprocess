rule doublet_scdblfinder:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample, "raw_feature_bc_matrix"),
        whitelist_path = lambda wc: os.path.join(RESULTS_DIR, "empty", wc.sample, "whitelist.txt"),
        dehash_path = lambda wc: os.path.join(RESULTS_DIR, "dehash", wc.sample, "barcode_metadata.csv")
    output:
        barcode_metadata = os.path.join(RESULTS_DIR, "doublet", "{sample}", "barcode_metadata.csv"),
        output = os.path.join(RESULTS_DIR, "doublet", "{sample}", "output.qs"),
        plot = os.path.join(RESULTS_DIR, "doublet", "{sample}", "plots.pdf"),
        session = os.path.join(RESULTS_DIR, "doublet", "{sample}", "session_info.txt"),
        whitelist = os.path.join(RESULTS_DIR, "doublet", "{sample}", "whitelist.txt")
    params:
        main = lambda wc: params_to_cli_args(config["doublet_scdblfinder"]),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, "doublet", wc.sample),
    log:
        os.path.join(RESULTS_DIR, "doublet", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --counts_dir {input.counts_dir} \
            --dehash_path {input.dehash_path} \
            --out_dir {params.out_dir} \
            --whitelist_path {input.whitelist_path} \
            {params.main} &> {log}
        """

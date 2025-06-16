rule qc_sc_sample:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample, "raw_feature_bc_matrix"),
        group_id_src_path = lambda wc: os.path.join(RESULTS_DIR, "dehash", wc.sample, "barcode_metadata.csv"),
        whitelist_path = lambda wc: os.path.join(RESULTS_DIR, "filter", wc.sample, "whitelist.txt")
    output:
        metrics = os.path.join(RESULTS_DIR, "qc_sc_sample", "{sample}", "metrics.csv"),
        session = os.path.join(RESULTS_DIR, "qc_sc_sample", "{sample}", "plots.pdf")
    params:
        main = lambda wc: params_to_cli_args(config["qc_sc_sample"]),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, "qc_sc_sample", wc.sample)
    log:
        os.path.join(RESULTS_DIR, "qc_sc_sample", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --counts_dir {input.counts_dir} \
            --group_id_src_path {input.group_id_src_path} \
            --whitelist_path {input.whitelist_path} \
            --out_dir {params.out_dir} \
            {params.main} &> {log}
        """

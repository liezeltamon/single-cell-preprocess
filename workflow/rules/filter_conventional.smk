rule filter_conventional:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample, DATA_DIR_SUFFIX),
        sample_id_src_path = lambda wc: os.path.join(RESULTS_DIR, "dehash", wc.sample, "barcode_metadata.csv"),
        whitelists = lambda wc: [
            os.path.join(RESULTS_DIR, "dehash", wc.sample, "whitelist.txt"),
            os.path.join(RESULTS_DIR, "doublet", wc.sample, "whitelist.txt")
        ]
    output:
        barcode_metadata = os.path.join(RESULTS_DIR, "filter", "{sample}", "barcode_metadata.csv"),
        plot = os.path.join(RESULTS_DIR, "filter", "{sample}", "plots.pdf"),
        session = os.path.join(RESULTS_DIR, "filter", "{sample}", "session_info.txt"),
        whitelist = os.path.join(RESULTS_DIR, "filter", "{sample}", "whitelist.txt")
    params:
        main = lambda wc: params_to_cli_args(config["filter_conventional"]),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, "filter", wc.sample)
    log:
        os.path.join(RESULTS_DIR, "filter", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --counts_dir {input.counts_dir} \
            --out_dir {params.out_dir} \
            --sample_id_src_path {input.sample_id_src_path} \
            --whitelists {input.whitelists} \
            {params.main} &> {log}
        """

rule empty_emptydrops:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample, "raw_feature_bc_matrix")
    output:
        blacklist = os.path.join(RESULTS_DIR, "empty", "{sample}", "blacklist.txt"),
        output = os.path.join(RESULTS_DIR, "empty", "{sample}", "output.qs"),
        plot = os.path.join(RESULTS_DIR, "empty", "{sample}", "plots.pdf"),
        session = os.path.join(RESULTS_DIR, "empty", "{sample}", "session_info.txt"),
        whitelist = os.path.join(RESULTS_DIR, "empty", "{sample}", "whitelist.txt")
    params:
        main = lambda wc: params_to_cli_args(config["empty_emptydrops"]),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, "empty", wc.sample),
    log:
        os.path.join(RESULTS_DIR, "empty", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --counts_dir {input.counts_dir} \
            --out_dir {params.out_dir} \
            {params.main} &> {log}
        """

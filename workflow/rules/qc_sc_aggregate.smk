rule qc_sc_aggregate:
    input:
        sample_metrics = expand(
            os.path.join(RESULTS_DIR, "qc_sc_sample", "{sample}", "metrics.csv"),
            sample=SAMPLES
        ),
        sample_plots = expand(
            os.path.join(RESULTS_DIR, "qc_sc_sample", "{sample}", "plots.pdf"),
            sample=SAMPLES
        )
    output:
        heatmaps = os.path.join(RESULTS_DIR, "qc_sc_aggregate", "heatmaps.pdf"),
        metrics = os.path.join(RESULTS_DIR, "qc_sc_aggregate", "metrics.csv")
    params:
        main = lambda wc: params_to_cli_args(config["qc_sc_aggregate"]),
        src_dir = os.path.join(RESULTS_DIR, "qc_sc_sample"),
        out_dir = os.path.join(RESULTS_DIR, "qc_sc_aggregate")
    log:
        os.path.join(RESULTS_DIR, "qc_sc_aggregate", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --src_dir {input.src_dir} \
            --out_dir {params.out_dir} \
            {params.main} &> {log}
        """

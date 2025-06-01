import os
# sbatch -p long --mem=50G --cpus-per-task=3 --output=log.out --error=log.err --wrap="snakemake -s rules/empty.smk --configfile envs/config.yml --cores 2 -j 2" &

# Parameters

PROJECT_DIR = os.path.abspath("../../../..")
DATA_DIR = os.path.join(PROJECT_DIR, "data", "cellranger")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "process_droplets_pipeline")

configfile: os.path.join("envs", "config.yml")

# Functions

def params_to_cli_args(param_dict):
    args = []
    for k, v in param_dict.items():
        flag = f"--{k}"
        if isinstance(v, bool):
            if v:
                args.append(flag)
        else:
            args.append(f"{flag} {v}")
    return " ".join(args)

def params_to_string(params_dict, parent=None):
    parts = []
    for k in sorted(params_dict):
        v = params_dict[k]
        if isinstance(v, float):
            v_str = f"{v:.4f}".replace('.', '')
        else:
            v_str = str(v)
        parts.append(f"{k}{v_str}")
    
    param_dirname = "_".join(parts)
    
    if parent:
        return f"{parent}/{param_dirname}"
    else:
        return param_dirname

PARAM_MAP = {
    name: params_to_string(params)
    for name, params in config.items()
}
PARAM_NAMES = list(PARAM_MAP.keys())
PARAM_DIRNAMES = list(PARAM_MAP.values())

SAMPLES = [d for d in os.listdir(DATA_DIR) if os.path.isdir(os.path.join(DATA_DIR, d))]

# Rules

rule all:
  input:
      expand(
          os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "{filename}"),
          param_dirname=PARAM_DIRNAMES,
          sample=SAMPLES,
          filename=[
              "blacklist.txt",
              "whitelist.txt",
              "output.qs",
              "plots.pdf",
              "session_info.txt"
          ]
      )
rule empty_emptydrops:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample, "raw_feature_bc_matrix")
    output:
        blacklist = os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "blacklist.txt"),
        whitelist = os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "whitelist.txt"),
        output = os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "output.qs"),
        plot = os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "plots.pdf"),
        session = os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "session_info.txt")
    params:
        main = lambda wc: params_to_cli_args(
            config[[k for k, v in PARAM_MAP.items() if v == wc.param_dirname][0]]
        ),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, wc.param_dirname, wc.sample),
    log:
        os.path.join(RESULTS_DIR, "{param_dirname}", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/{rule}.R \
            --counts_dir {input.counts_dir} \
            --out_dir {params.out_dir} \
            {params.main} &> {log}
        """
        
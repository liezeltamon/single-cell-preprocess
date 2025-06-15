def params_to_cli_args(param_dict):
    args = []
    for k, v in param_dict.items():
        flag = f"--{k}"
        if isinstance(v, bool):
            if v:
                args.append(flag)
        elif isinstance(v, list):
            args.append(f"{flag} " + " ".join(map(str, v)))
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

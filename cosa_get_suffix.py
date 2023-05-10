def cosa_get_suffix(anaerobic: bool, expanded: bool, c_source: str="glucose", fixed_nadx_change: int=0) -> str:
    suffix = ""
    if anaerobic:
        suffix += "_anaerobic"
    else:
        suffix += "_aerobic"
    if expanded:
        suffix += "_expanded"

    if c_source != "glucose":
        suffix += f"_{c_source}"

    if fixed_nadx_change != 0:
        suffix += f"_nadz_change_{fixed_nadx_change}"

    return suffix

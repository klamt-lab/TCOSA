def cosa_get_suffix(anaerobic: bool, expanded: bool) -> str:
    suffix = ""
    if anaerobic:
        suffix += "_anaerobic"
    else:
        suffix += "_aerobic"
    if expanded:
        suffix += "_expanded"
    return suffix

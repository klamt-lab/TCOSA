def cosa_get_suffix(anaerobic: bool, expanded: bool, c_source: str="glucose") -> str:
    suffix = ""
    if anaerobic:
        suffix += "_anaerobic"
    else:
        suffix += "_aerobic"
    if expanded:
        suffix += "_expanded"

    if c_source != "glucose":
        suffix += f"_{c_source}"

    return suffix

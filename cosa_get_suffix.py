"""Defines function for adding a cosa results folder suffix."""

def cosa_get_suffix(anaerobic: bool, expanded: bool, c_source: str="glucose", fixed_nadx_change: int=0) -> str:
    """Adds the TCOSA results folder suffix appropriately.

    Args:
        anaerobic (bool): _description_
        expanded (bool): _description_
        c_source (str, optional): _description_. Defaults to "glucose".
        fixed_nadx_change (int, optional): _description_. Defaults to 0.

    Returns:
        str: _description_
    """
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

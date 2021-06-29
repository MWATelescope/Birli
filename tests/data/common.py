def chunk(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip(*args)

def get_gpufits_num_scans(hdu_len, corr_type):
    if corr_type == "MWAX":
        return (hdu_len - 1) // 2
    return hdu_len - 1
def get_copy_call(a):
    """
    Heler Function - Determine copy call for a value.
    """

    if a is None:
        return float('nan')

    if math.isnan(a):
        return float('nan')

    a_val = math.log2(float(a)+0.000001)
    if a_val < 0.5210507:
        return 'deep del'
    elif a_val < 0.7311832:
        return 'het loss'
    elif a_val < 1.214125:
        return 'diploid'
    elif a_val < 1.422233:
        return 'gain'
    else:
        return 'amp'

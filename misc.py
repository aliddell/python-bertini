def dps(s):
    #return max(DPS, len(s))
    s = str(s) if not isinstance(s, str) else s
    return len(s)


def striplines(lines, nonempty=True):
    if nonempty:
        return [l.strip() for l in lines if l != '\n']
    else:
        return [l.strip() for l in lines]
from __future__ import absolute_import, print_function


def striplines(lines, nonempty=True):
    if nonempty:
        return [l.strip() for l in lines if l != '\n']
    else:
        return [l.strip() for l in lines]
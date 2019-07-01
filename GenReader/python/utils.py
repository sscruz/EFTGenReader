import re

# Match strings using one or more regular expressions
def regex_match(lst,regex_lst,verbose=0):
    # NOTE: We don't escape any of the regex special characters!
    matches = []
    if len(regex_lst) == 0:
        return lst[:]
    if verbose:
        for p in regex_lst: print "rgx:",r"%s" % (p)
    for s in lst:
        for pat in regex_lst:
            m = re.search(r"%s" % (pat),s)
            if m is not None:
                matches.append(s)
                break
    return matches
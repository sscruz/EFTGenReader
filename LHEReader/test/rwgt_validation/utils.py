import os
import re

# Match strings using one or more regular expressions
def regex_match(lst,regex_lst):
    # NOTE: We don't escape any of the regex special characters!
    # TODO: Add whitelist/blacklist option switch
    matches = []
    if len(regex_lst) == 0:
        return lst[:]
    for s in lst:
        for pat in regex_lst:
            m = re.search(r"%s" % (pat),s)
            if m is not None:
                matches.append(s)
                break
    return matches
    
# Returns all info objects which matches the info_tag
# Note: Matching is done based on the tail of the path defined by the tag, so multiple results are possible
def getInfoByTag(lst,info_tag):
    sub_lst = []
    for info in lst:
        head,tail = os.path.split(info['tag'])
        if info_tag == tail:
            sub_lst.append(info)
    return sub_lst

# Get list of all directories which pass the whitelists
def getDirectories(path,p_wl=[],c_wl=[],r_wl=[]):
    dir_list = []
    for fd in os.listdir(path):
        if not os.path.isdir(os.path.join(path,fd)):
            continue
        arr = fd.split('_')
        if len(arr) != 4:
            print "[WARNING] Bad name: %s" % (fd)
            continue
        p,c,r = arr[1],arr[2],arr[3]
        if len(p_wl) > 0 and not p in p_wl:
            continue
        elif len(c_wl) > 0 and not c in c_wl:
            continue
        elif len(r_wl) > 0 and not r in r_wl:
            continue
        dir_list.append(os.path.join(path,fd))
    return dir_list

# Group the directories based on process
def groupByProcess(path,tag,grp='',p_wl=[],c_wl=[],r_wl=[]):
    grp_dirs = {}   # {(tag,ttH,grp): [path1,path2,...]}
    for fd in os.listdir(path):
        if not os.path.isdir(os.path.join(path,fd)):
            continue
        arr = fd.split('_')
        if len(arr) != 4:
            print "[WARNING] Bad name: %s" % (fd)
            continue
        p,c,r = arr[1],arr[2],arr[3]
        if len(p_wl) > 0 and not p in p_wl:
            continue
        elif len(c_wl) > 0 and not c in c_wl:
            continue
        elif len(r_wl) > 0 and not r in r_wl:
            continue
        #if len(regex_match([p],p_wl)) == 0:
        #    continue
        #elif len(regex_match([c],c_wl)) == 0:
        #    continue
        #elif len(regex_match([r],r_wl)) == 0:
        #    continue
        key = (tag,p,grp)
        if not grp_dirs.has_key(key):
            grp_dirs[key] = []
        grp_dirs[key].append(os.path.join(path,fd))
    return grp_dirs

# Group the directories based on process and coeff tags
def groupByCoefficient(path,tag,p_wl=[],c_wl=[],r_wl=[]):
    grp_dirs = {}   # {(tag,ttH,grp): [path1,path2,...]}
    for fd in os.listdir(path):
        if not os.path.isdir(os.path.join(path,fd)):
            continue
        arr = fd.split('_')
        if len(arr) != 4:
            print "[WARNING] Bad name: %s" % (fd)
            continue
        p,c,r = arr[1],arr[2],arr[3]
        if len(p_wl) > 0 and not p in p_wl:
            continue
        elif len(c_wl) > 0 and not c in c_wl:
            continue
        elif len(r_wl) > 0 and not r in r_wl:
            continue
        #key = "%s_%s" % (p,c)
        key = (tag,p,c)
        if not grp_dirs.has_key(key):
            grp_dirs[key] = []
        grp_dirs[key].append(os.path.join(path,fd))
    return grp_dirs
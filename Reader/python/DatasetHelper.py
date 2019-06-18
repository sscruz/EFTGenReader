import os
import subprocess
import json

class DatasetHelper(object):
    def __init__(self,fn=None):
        self.datasets = {}
        self.root_redirect = "root://ndcms.crc.nd.edu/"
        self.hadoop_protocol = "file://"
        if fn: self.load(fn)

    # Load the dataset info from a json file
    def load(self,fn,update=False):
        if not os.path.exists(fn):
            return
        if not update:
            self.datasets = {}
        with open(fn) as f:
            d = json.load(f)
            for k,ds in d.iteritems():
                if not self.datasets.has_key(k):
                    self.datasets[k] = DSContainer()
                self.datasets[k].setData(**ds)

    def getFiles(self,name):
        if not self.datasets.has_key(name):
            return []
        ds = self.datasets[name]
        if ds.getData('on_das'):
            lst = self.findOnDas(ds)
        else:
            lst = self.findOnHadoop(ds)
        return lst

    def getData(self,name,data_field):
        if not self.datasets.has_key(name):
            return None
        return self.datasets[name].getData(data_field)

    def findOnDAS(self,ds):
        lst = []
        das_query = 'file dataset=%s | grep file.name, file.size, file.nevents' % (ds.getData('dataset'))
        ret = subprocess.check_output(['dasgoclient','--query',query])
        ret = ret.split('\n')
        max_files = 10
        for idx,l in enumerate(ret):
            if idx >= max_files: break
            fn,sz,evts = l.split()
            lst.append(self.root_redirect + fn)
        return lst

    def findOnHadoop(self,ds):
        lst = []
        loc = ds.getData('loc')
        if not os.path.exists(loc):
            print "Unknown fpath: %s" % (loc)
            return []
        for fn in os.listdir(loc):
            if not ".root" in fname: continue
            fpath = self.hadoop_protocol + os.path.join(loc,fn)
            lst.append(fpath)

class DSContainer(object):
    def __init__(self,**kwargs):
        self.data = {
            'dataset': '',
            'loc': None,
            'on_das': False,
            'is_eft': False,
            'central_xsec': 0.0
        }

    def setData(self,**kwargs):
        for k,v in kwargs.iteritems():
            if not self.data.has_key(k):
                continue
            self.data[k] = v

    def getData(self,data_field):
        if not self.data.has_key(data_field):
            return None
        return self.data[data_field]
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
        # update: If set to true, will simply overwrite existing datasets and add new ones
        if not os.path.exists(fn):
            return
        if not update:
            self.datasets = {}
        with open(fn) as f:
            d = json.load(f)
            for k,ds in d.iteritems():
                self.updateDataset(k,**ds)

    # Save the dataset info to a json file
    def save(self,fn):
        with open(fn,'w') as f:
            d = self.toDict()
            json.dump(d,f,indent=2)

    # Remove all datasets
    def clear(self):
        self.datasets.clear()

    # Create/Update/Modify a dataset
    def updateDataset(self,name,**kwargs):
        if not self.datasets.has_key(name):
            self.datasets[name] = DSContainer()
        self.datasets[name].setData(**kwargs)

    # Remove a specific dataset
    def removeDataset(self,name):
        if self.datasets.has_key(name):
            del self.datasets[name]

    # Returns a list of root file locations corresponding to the specified dataset
    def getFiles(self,name):
        if not self.datasets.has_key(name):
            return []
        ds = self.datasets[name]
        if ds.getData('on_das'):
            lst = self.findOnDAS(ds)
        else:
            lst = self.findOnHadoop(ds)
        return lst

    # Passthrough to underlying dataset getter
    def getData(self,name,data_field):
        if not self.datasets.has_key(name):
            return None
        return self.datasets[name].getData(data_field)

    # Try and find dataset root files from DAS
    def findOnDAS(self,ds):
        lst = []
        das_query = 'file dataset=%s | grep file.name, file.size, file.nevents' % (ds.getData('dataset'))
        ret = subprocess.check_output(['dasgoclient','--query',das_query])
        ret = ret.split('\n')
        max_files = 10      # currently arbitrary
        for idx,l in enumerate(ret):
            if idx >= max_files: break
            fn,sz,evts = l.split()
            lst.append(self.root_redirect + fn)
        return lst

    # Try and find dataset root files on the local hadoop cluster
    def findOnHadoop(self,ds):
        lst = []
        loc = str(ds.getData('loc'))    # converts from a unicode string object to a normal str
        if not os.path.exists(loc):
            print "Unknown fpath: %s" % (loc)
            return []
        idx = 0
        max_files = 100     # currently arbitrary
        for fn in os.listdir(loc):
            if not ".root" in fn: continue
            if idx >= max_files: break
            fpath = self.hadoop_protocol + os.path.join(loc,fn)
            lst.append(fpath)
            idx += 1
        return lst

    # Returns all datasets in a format suitable to be passed to json.dump()
    # NOTE: The nested dictionaries will be references to the ones stored in the actual
    #       DSContainer() objects
    def toDict(self):
        d = {}
        for k in self.datasets.keys():
            d[k] = self.datasets[k].toDict()
        return d

class DSContainer(object):
    def __init__(self,**kwargs):
        self.__data = {
            'dataset': '',
            'loc': None,
            'on_das': False,
            'is_eft': False,
            'central_xsec': 0.0
        }
        self.setData(**kwargs)

    # Generic setter
    def setData(self,**kwargs):
        for k,v in kwargs.iteritems():
            if not self.data.has_key(k):
                continue
            self.__data[k] = v

    # Generic getter
    def getData(self,data_field):
        if not self.data.has_key(data_field):
            return None
        return self.__data[data_field]

    # Convert self to dictionary suitable to be passed to json.dump()
    def toDict(self):
        return self.__data
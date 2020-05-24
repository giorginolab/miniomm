import configparser
from collections import OrderedDict

root_section = "ROOT"


class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            # super(MultiOrderedDict, self).__setitem__(key, value)
            super().__setitem__(key, value)


class Config():
    def __init__(self, fn):
        self.cp = configparser.RawConfigParser(dict_type=MultiOrderedDict,
                                               delimiters=[' ', "\t", '='],
                                               strict=False)
        with open(fn) as f:
            self.cp.read_string('['+root_section+"]\n"+f.read())

    def get(self, k):
        return self.cp.get(root_section, k)

    def __getattr__(self, k):
        return self.get(k)

    def __contains__(self, k):
        return k in self.cp[root_section]

    def getboolean(self, k):
        v = self.__getattr__(k).lower()
        return (v == "on" or v == "true")


if __name__ == "__main__":
    import sys

    c = Config("scratch/pdir_amber/input")
    c.cp.write(sys.stdout)

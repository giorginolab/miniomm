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
        self.usedKeys = {}
        self.cp = configparser.RawConfigParser(dict_type=MultiOrderedDict,
                                               delimiters=[' ', "\t", '='],
                                               strict=False)
        with open(fn) as f:
            self.cp.read_string('['+root_section+"]\n"+f.read())

    def get(self, k):
        self.usedKeys[k.lower()] = 1
        return self.cp.get(root_section, k)

    def __getattr__(self, k):
        return self.get(k)

    def __contains__(self, k):
        self.usedKeys[k.lower()] = 1
        return k in self.cp[root_section]

    def getboolean(self, k):
        self.usedKeys[k.lower()] = 1
        v = self.__getattr__(k).lower()
        if v == "on" or v == "true":
            return True
        elif v == "off" or v == "false":
            return False
        else:
            raise ValueError(f"Configuration key {k} not a boolean")

    def unusedKeys(self):
        kl1 = self.cp.items(root_section)
        kl2 = [x[0].lower() for x in kl1]
        kl3 = set(kl2)
        sdiff = kl3.difference(self.usedKeys.keys())
        return sdiff


if __name__ == "__main__":
    import sys

    c = Config("scratch/pdir_amber/input")
    c.cp.write(sys.stdout)

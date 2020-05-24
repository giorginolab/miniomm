import configparser


class Config():
    def __init__(self,fn):
        self.cp = configparser.ConfigParser(delimiters=[' ',"\t",'='])
        with open(fn) as f:
            self.cp.read_string("[root]\n"+f.read())


if __name__ == "__main__":
    import sys
    
    c=Config("scratch/pdir_amber/input")
    c.cp.write(sys.stdout)
    


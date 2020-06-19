import simtk.openmm as mm

from miniomm.miniomm import run_omm
import miniomm.util 
import miniomm


def _banner():
    return f"""
            _         _                              
 _ __ ___  (_) _ __  (_)  ___   _ __ ___   _ __ ___  
| '_ ` _ \ | || '_ \ | | / _ \ | '_ ` _ \ | '_ ` _ \ 
| | | | | || || | | || || (_) || | | | | || | | | | |
|_| |_| |_||_||_| |_||_| \___/ |_| |_| |_||_| |_| |_|
                                               {miniomm.__version__}           
A minimalistic OpenMM MD frontend.   
https://github.com/giorginolab/miniomm
"""


def main():
    from optparse import OptionParser
    parser = OptionParser()
    platformNames = [mm.Platform.getPlatform(
        i).getName() for i in range(mm.Platform.getNumPlatforms())]
    parser.add_option('--input', dest='input',
                      default="input", help='name of the input file')
    parser.add_option('--platform', dest='platform',
                      choices=platformNames, help='name of the platform to benchmark')
    parser.add_option('--device', default=None, dest='device',
                      help='device index for CUDA or OpenCL')
    parser.add_option('--precision', dest='precision', choices=('single', 'mixed', 'double'),
                      help='precision mode for CUDA or OpenCL: single, mixed, or double')
    parser.add_option('--hours', default='11.5', dest='run_hours', type='float',
                      help='target simulation length in hours [default: 11.5]')

    (options, args) = parser.parse_args()

    if len(args) > 0:
        print("Remaining args: "+" ".join(args))

    print(_banner())
    run_omm(options)


if __name__ == "__main__":
    main()
    

import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "miniomm",
    version = "0.0.10",
    author = "Toni Giorgino",
    author_email = "toni.giorgino@gmail.com",
    description = ("A simple OpenMM wrapper for common-case MD runs"),
    license = "LGPL",
    keywords = "molecular dynamics",
    url = "https://github.com/giorginolab/miniomm",
    packages=['miniomm', 'tests'],
    entry_points = {
        'console_scripts': ['miniomm=miniomm.miniomm:main']
    },
    scripts = ['bin/sbatch_many'],
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)"
    ],
    # install_requires = [ 'simtk.openmm', 'simtk.unit' ],
)

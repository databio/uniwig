from setuptools import setup, Extension
import os
import sys

PACKAGE = "smugwig"

# Additional keyword arguments for setup().
extra = {}

# It might at some point be required to change the compiler to g++. Here's how:
# os.environ["CC"] = "g++"

# Ordinary dependencies
DEPENDENCIES = []
with open("requirements/requirements-all.txt", "r") as reqs_file:
    for line in reqs_file:
        if not line.strip():
            continue
        #DEPENDENCIES.append(line.split("=")[0].rstrip("<>"))
        DEPENDENCIES.append(line)

if sys.version_info >= (3, ):
    extra["use_2to3"] = True
extra["install_requires"] = DEPENDENCIES


with open("{}/_version.py".format(PACKAGE), 'r') as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")

# Handle the pypi README formatting.
try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
except(IOError, ImportError, OSError):
    long_description = open('README.md').read()

setup(
	name = PACKAGE,
	packages = [PACKAGE],
	version = version,
	description="A fast converter for bam files to wiggle",
    long_description=long_description,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],    
	maintainer = "Nathan Sheffield",
	author = u"Nathan Sheffield, Jason Smith",
	license="BSD2",
	ext_modules = [Extension('smugwigc', sources=['src/smugwig.cpp'])],
    python_requires='>=3.4',
    entry_points = {
    "console_scripts": ['smugwigi = smugwig.swinternal:main', 
    					'smugwig = smugwig.smugwig:main'],
	},
	**extra
)

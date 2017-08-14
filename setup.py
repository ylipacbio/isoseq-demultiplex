import os
import re
from setuptools import find_packages

try:
    from setuptools import setup
    from distutils.core import setup
except ImportError:
    pass

version = __import__('debarcode').get_version()

_REQUIREMENTS_FILE = 'REQUIREMENTS.txt'
_README = 'README.md'


def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


def _get_description(file_name):
    with open(file_name, 'r') as f:
        _long_description = f.read()
    return _long_description


def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    rx = re.compile('^[A-z]')
    requirements = [l for l in lines if rx.match(l) is not None]
    if "READTHEDOCS" in os.environ:
        requirements = [r for r in requirements if not "pbcore" in r]
    return requirements

setup(
    name='isoseq-demultiplex',
    version=version,
    package_dir={'debarcode': 'debarcode'},
    packages=find_packages('.'),
    license='BSD',
    author='yli',
    author_email='yli@pacificbiosciences.com',
    description='Trace zmws of barcoded isoseq sample to primers based on clustering',
    setup_requires=['nose>=1.0'],
    # Maybe the pbtools-* should really be done in a subparser style
    entry_points={'console_scripts': [
        'zmw-to-consensus-primer = debarcode.zmw_to_consensus_primer:main',
        'cluster-to-consensus-primer = debarcode.cluster_to_consensus_primer:main'
    ]},
    install_requires=_get_requirements(_get_local_file(_REQUIREMENTS_FILE)),
    tests_require=['nose'],
    long_description=_get_description(_get_local_file(_README)),
    classifiers=['Development Status :: 4 - Beta'],
    include_package_data=True,
    zip_safe=False
)

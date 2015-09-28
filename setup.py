#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup, Command
except ImportError:
    from distutils.core import setup, Command

# Run pre-built py.tests as described at
# https://pytest.org/latest/goodpractises.html#integrating-with-distutils-python-setup-py-test
class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

with open('DESCRIPTION.rst') as f:
    long_description = f.read()

requirements = [
    'requests>=2.7.0',
    'progressbar>=2.3',
    'astropy>=1.0.1',
    'fitsio>=0.9.7',
    'numpy>=1.9.2',
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='bossdata',
    version='0.2.7',
    description='Tools to access SDSS spectroscopic data.',
    long_description=long_description,
    author='bossdata developers',
    author_email='dkirkby@uci.edu',
    url='https://github.com/dkirkby/bossdata',
    packages=[
        'bossdata',
    ],
    package_dir={'bossdata':
                 'bossdata'},
    scripts = [
        'bin/bossfetch',
        'bin/bossquery',
        'bin/bossplot',
    ],
    #include_package_data=True,
    #zip_safe=False,
    install_requires=requirements,
    license='MIT',
    keywords='bossdata',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    cmdclass = {'test': PyTest},
    #test_suite='tests',
    #tests_require=test_requirements
)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('DESCRIPTION.rst') as f:
    long_description = f.read()

requirements = [
    'requests',
    'progressbar2',
    'numpy',
    'astropy',
    'fitsio',
    'pydl',
    'six'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='bossdata',
    version='0.3.1',
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
        'bin/bossraw',
        'bin/bosssky',
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
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest']
)

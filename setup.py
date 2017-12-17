#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('DESCRIPTION.rst') as f:
    long_description = f.read()

requirements = [
    'requests>=2.14',
    'progressbar2>=3.34',
    'numpy>=1.13',
    'astropy>=2.0.1',
    'fitsio>=0.9.11',
    'pydl>=0.6.0',
    'six>=1.10'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='bossdata',
    version='0.3.0',
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
    setup_requires=['pytest-runner'],
    tests_require=['pytest']
)

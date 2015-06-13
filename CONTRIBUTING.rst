============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/dkirkby/bossdata/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to fix it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

bossdata could always use more documentation, whether as part of the
official bossdata docs, in docstrings, or even on the web in blog posts,
articles, and such.  We use the `sphinx napolean extension <http://sphinx-doc.org/latest/ext/napoleon.html>`_ and write google-style docstrings. Some helpful tips:

* Use ```text <http://url>`_`` to embed external links (don't forget the space!)
* Add ``.. _scriptname:`` before the heading for new scripts in `bin/scripts.rst`.  You can refer to these from other markup as ``:ref:`scriptname```.
* Refer to another markup document `otherdoc.rst` as ``:doc:`otherdoc```.
* Add cross references to locally defined API entities using:

 * classes ``:class:`bossdata.module.Class```
 * methods ``:meth:`bosdata.module.Class.method```
 * functions ``:func:`bossdata.module.func```

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/dkirkby/bossdata/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `bossdata` for local development.

1. Fork the `bossdata` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/bossdata.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv bossdata
    $ cd bossdata/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox::

    flake8 --doctests --exclude bossdata/bits.py --max-line-length 95 bossdata
    py.test --doctest-modules --verbose bossdata

   Note that `--doctest-modules` will require that all external modules imported from our modules are installed, so omit that option if you only want to run the unit tests.

   To get flake8 and tox, just pip install them into your virtualenv.  If you don't already have a python-2.6 environment installed, try this if you are using conda::

    conda create -n py26 python=2.6 anaconda

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Version Update Checklist
------------------------

1. Start a new release candidate branch, e.g., 0.2.1rc
2. Update the `version` in `setup.py`
3. Update the `__version__` in `__init__.py`
4. Create a pull request on github for this branch (its ok that it doesn't have any new code yet).
5. Iterate on changes.
6. Add a brief description of the changes to `HISTORY.rst`
7. Push changes to github, which will trigger a Travis integration test of the release-candidate branch.
8. Merge the pull request.
9. Update local master and tag the new version::

    git fetch
    git checkout master
    git pull
    git tag 0.2.1
    git push --tags
    git branch -d 0.2.1rc

9. Submit the changes to pypi::

    python setup.py sdist bdist_wheel upload

New External Depencency Checklist
------------------------

These steps are not required for modules that are included with the python standard library.

1. Add to `MOCK_MODULES` in `docs/conf.py`.
2. Add the actual version being used to `requirements.txt`
3. Add to the `requirements` list in `setup.py`
4. Mention in `docs/installation.rst`

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.6 and 2.7. Check
   https://travis-ci.org/dkirkby/bossdata/pull_requests
   and make sure that the tests pass for all supported Python versions.

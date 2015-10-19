============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Issues
~~~~~~~~~~~~~

Report issues on our `issues page <https://github.com/dkirkby/bossdata/issues>`_. First check that if your issue is already addressed.  If so, feel free to join its conversation and add any relevant information from your experience.  If this is a new issue, click the ``New Issue`` button to describe it, including:

* The type of data you are trying to access (BOSS, SEQUELS, ...)
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Propose a New Feature
~~~~~~~~~~~~~~~~~~~~~

You can also open a new issue to propose a new feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Work on Issues
~~~~~~~~~~~~~~

Look through the `open issues <https://github.com/dkirkby/bossdata/issues>`_ for areas where we currently need help from developers like you. If you find an issue that you are willing to contribute to, start by joining its conversation and tell us about your ideas.

Write Documentation
~~~~~~~~~~~~~~~~~~~

bossdata could always use more documentation, whether as part of the
official bossdata docs, in docstrings, or even on the web in blog posts,
articles, and such.

We use the `sphinx napolean extension <http://sphinx-doc.org/latest/ext/napoleon.html>`_ and write google-style docstrings. Some helpful tips:

* Use ```text <http://url>`_`` to embed external links (don't forget the space!)
* Add ``.. _scriptname:`` before the heading for new scripts in `bin/scripts.rst`.  You can refer to these from other markup as ``:ref:`scriptname```.
* Refer to another markup document `docs/otherdoc.rst` as ``:doc:`/otherdoc```.
* Add cross references to locally defined API entities using:

 * classes ``:class:`bossdata.module.Class```
 * methods ``:meth:`bosdata.module.Class.method```
 * functions ``:func:`bossdata.module.func```

* You can override the default link text by changing ``:role:`target``` to ``:role:`text <target>```.

Get Started!
------------

Ready to contribute? Here's how to set up ``bossdata`` for local development.

1. Fork the ``bossdata`` repo on GitHub.
2. Clone your fork locally::

    git clone git@github.com:your_name_here/bossdata.git

3. Install your local copy for local development::

    cd bossdata/
    python setup.py develop --user

   To later revert back to a system-installed version of the package, un-install your development install using::

    python setup.py develop --user --uninstall

4. Create a branch for local development::

    git checkout -b '#nnn'
    git push -u origin '#nnn'

   where ``nnn`` is the number of the issue you are working on (quotes are required because of the ``#`` symbol in the branch name). Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the unit tests::

    flake8 --doctests --exclude bossdata/bits.py --max-line-length 95 bossdata
    py.test --doctest-modules --verbose bossdata

   Note that `--doctest-modules` will require that all external modules imported from our modules are installed, so omit that option if you only want to run the unit tests.  If you don't already have flake8, you can pip install it.

6. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin '#nnn'

7. Submit a `pull request <https://github.com/dkirkby/bossdata/pulls>`_.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests, if appropriate.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in HISTORY.rst.
3. The pull request should work for Python 2.6 and 2.7. Check
   https://travis-ci.org/dkirkby/bossdata/pull_requests
   and make sure that the tests pass for all supported Python versions.

Version Update Checklist
------------------------

#. Start a new release candidate branch, e.g::

    git checkout -b 0.2.1rc
    git push -u origin 0.2.1rc

#. Update the ``version`` in ``setup.py``
#. Update the ``__version__`` in ``__init__.py``
#. Add a brief description of the changes to ``HISTORY.rst`` and update the ``What's New`` section of ``DESCRIPTION.rst`` (which is what pypi will display for this release). You can get a list of merges to master since the last tagged release using::

    git log --oneline --merges `git describe --tags --abbrev=0`..HEAD

#. Push changes to github, which will trigger a Travis integration test of the release-candidate branch.
#. Create a pull request on github for this branch and ask someone else to review it and give feedback.
#. Merge the pull request.
#. Update local master and tag the new version, e.g::

    git fetch
    git checkout master
    git pull
    git tag 0.2.1
    git push --tags
    git branch -d 0.2.1rc

#. Submit the changes to pypi::

    python setup.py sdist bdist_wheel upload

#. Update the ``version`` in ``setup.py`` and ``__version__`` in ``__init__.py`` to indicate that master is under development, e.g. to ``0.2.2dev``.
#. Reset the ``What's New`` section of ``DESCRIPTION.rst`` and add a new entry at the bottom of ``HISTORY.rst``, e.g::

    0.2.2 (unreleased)
    ------------------

    * No changes yet.

#. Update master so that new topic branches will include these changes, e.g::

    git add setup.py bossdata/__init__.py HISTORY.rst DESCRIPTION.rst
    git commit -m 'Start development on version 0.2.2'
    git push

New External Depencency Checklist
---------------------------------

These steps are not required for modules that are included with the python standard library.

1. Add to `MOCK_MODULES` in `docs/conf.py`.
2. Add the actual version being used to `requirements.txt`
3. Add to the `requirements` list in `setup.py`
4. Mention in `docs/installation.rst`

.. _code-development:

Code development
================

This page includes some useful information and pointers for people willing
to contribute changes to ixpeobssim. We assume that you have read
the :ref:`installation` page before landing here.


.. warning::
   This page is fairly Linux-centric, in that most of the development, this
   days, happens under GNU-Linux. That said, all the components and tools that
   we use are instrinsically cross-platform, and ixpeobssim is known to be
   working on Windows and Mac. Feel free to edit this page making it more
   friendly for Windows and Mac users.


Up and running with github
--------------------------

`git <http://git-scm.com/>`_ is a distributed version control system and
`github <https://github.com>`_ is the web hosting service that we use to develop
the public version of `ixpeobssim <https://github.com/lucabaldini/ixpeobssim>`_.
`Here <http://git-scm.com/doc>`_ is the entry point for the git documentation,
in case you want to have a feeling of what git is doing and how
you use it.

Mind that, in order to be able to push back changes to the remote repository
you will need to tell git on your machine who you are, i.e.:

.. code-block:: bash

    git config --global user.name "Your Name"
    git config --global user.email you@example.com


Cloning the repository
----------------------

In order to clone the repository, go to the webpage
`ixpeobssim <https://github.com/lucabaldini/ixpeobssim>`_, and click on the top
right "fork" icon. Create your fork, then clone it on your local device, by
typing:

.. code-block:: bash

    git clone git@github.com:github_username/ixpeobssim.git

(mind this will create an ixpeobssim folder in the current directory, so cd to
the appropriate place in your filesystem first).

If you get an error message along the lines of

.. code-block:: bash

   Permission denied (publickey).
   fatal: Could not read from remote repository.
   Please make sure you have the correct access rights and the repository exists.

that simply means that you have to exchange your public SSH key with the server.
In order to do that, click on your github profile icon on the top-right of
the github webpage, select "settings", "SSH and GPG keys", "New SSH key" (top right)
and paste in the form the content of the local (i.e. on
the machine you are cloning the repository into) ~/.ssh/id_rsa.pub file.

If you don’t have a public ssh key, you can generate it by typing

.. code-block:: bash

   ssh-keygen

(press ENTER a couple of times and here is your public key in ~/.ssh/id_rsa.pub)


Basic git workflow
------------------

The ``ixpeobssim`` public repository is intentionally protected,
meaning that nobody is allowed to push changes directly to it.

Everybody can merge changes onto the public repository,
via pull requests provided that there's a least one approval.
This scheme makes it for a fairly horizontal development approach, where
everybody can contribute changes more or less independently, but forces
people to do so in a coordinate fashion, and gives eveybody else a chance to
look at the code before changes are actually merged.

The basic workflow we want to stick to is essentially the following. Whenever
you are ready to start making a set of modifications, click on the "contribute"
icon on your fork webpage (i.e. https://github.com/github_username/ixpeobssim),
(on top of the "code" list), then "open pull request"; fill in the request and open it.

It is recommended to open a pull request from a branch of your fork,
rather than from the "main", in order to be able to work on more parallel tasks.
Create a new branch to work into and check it out (if you haven't already done so):

.. code-block:: bash

   git branch fixing_something
   git checkout fixing_something

It goes without saying that it is highly recommended to name the branch
making clear its intent (e.g., mybranch is not a very expressive name).

At this point you are in the new branch, and you can start doing your
modifications. Make sure your modifications do not break existing unit tests
(scroll down below for more information about that) and, if you are writing
brand new code, consider adding more unit tests covering the new territory.
Once you're happy with the changes, commit them

.. code-block:: bash

   git commit -m "Some expressive message" file1 file2 ... filen
   git push

Mind that the first time you push on the new branch you will get an
error message along the lines of

.. code-block:: bash

   git push

      fatal: The current branch fixing_something has no upstream branch.
      To push the current branch and set the remote as upstream, use
      git push --set-upstream origin fixing_something

Follow the instructions and you should be all set.

Once you are done with your consistent set of modifications, go ahead on the
repository web interface and create a pull request.
Click on the menu icon top left of the code list, on your fork webpage,
in order to select the right branch you want to make a pull request from
(default in this menu is "main"); then create and open your pull request,
as described above.
You're all set! Wait for the comments of the reviewer, and finally merge
the branch on the master (or, even better, have somebody else doing it for you).

.. _coding-guidelines:

Coding guidelines
-----------------

Though we'll never be able to follow any set of coding conventions religiously,
`PEP 0008 <https://www.python.org/dev/peps/pep-0008/>`_ is our starting point.
Take a second to give a look to this short recap of the most salient guidelines:

* Use 4 spaces for indentation level (no TABS).
* Limit all lines to 79 characters.
* Surround top-level function and class definitions with two blank lines.
  Method definitions inside a class are surrounded by a single blank line.
  Use blank lines in functions, sparingly, to indicate logical sections.
* Use one import per line, right at the top of the module.
* Use single quote characters for strings and double quotes characters for
  triple-quoted strings.
* Avoid extraneous white spaces, and especially avoid more than one space
  around an assigment.
* Don't use spaces around the `=` sign when used to indicate a keyword argument
  or a default parameter value.
* Modules should have short, all-lowercase names.
* Class names should normally use the CapWords convention (for ixpeobssim
  starting with a `x`).
* Function and member names should be lowercase, with words separated by
  underscores as necessary to improve readability.
* Constants are usually defined on a module level and written in all capital
  letters with underscores separating words.
* Always use a `def` statement instead of an assignment statement that binds a
  `lambda` expression directly to an identifier.

An example module, illustrating the basic guidelines, is available on the
repository at :repourl:`ixpeobssim/utils/codestyle.py`.


Documenting the code
--------------------

We use `sphinx <http://sphinx-doc.org/#>`_ to generate the ixpeobssim
documentation (which is the same big projects like Scipy, astropy and Python
itself are using). We use the `Napoleon
<https://sphinxcontrib-napoleon.readthedocs.org/en/latest/>`_ extension in the
Numpy flavor, and creating inline documentation essentially boils down to

* providing suitable docstrings with the appropriate syntax in the source files;
* creating suitable .rst files in the `doc/modules` folder.


In addition to `Napoleon`, you also will need `programoutput` and
`sphinx_rtd_theme` sphinx extensions. You can easily get them with ``pip`` running:

.. code-block:: bash

    python -m pip install sphinxcontrib-napoleon
    python -m pip install sphinxcontrib-programoutput
    python -m pip install sphinx_rtd_theme

Make sure also to have on your machine the `dvipng
<http://savannah.nongnu.org/projects/dvipng/>`_ package able to render math
equations via ``LaTeX``.

It won't take more than a few minutes to get aquainted to the basic rules,
and the :repourl:`ixpeobssim/utils/codestyle.py` module, along with its fellow
:repourl:`doc/modules/utils.codestyle.rst` file, provide a minimal working
example that, compiled with sphinx, would be rendered like
:mod:`ixpeobssim.utils.codestyle`.

You can compile and view the ixpeobssim documentation locally by doing

.. code-block:: bash

    cd docs
    make htmlall
    htmlview _build/html/index.html

which is useful to make sure everything is in order when writing and
documenting code.

Documentation is available online:
`<https://ixpeobssim.readthedocs.io/en/latest/overview.html>`

.. warning::
   We should update this section once the documentation is uploaded on the
   wiki and we have made up our mind about the access details.


Unit testing
------------

We use the Python `unittest <https://docs.python.org/2/library/unittest.html>`_
module for the purpose (the documentation includes a whole bunch of good
examples). While, again, we'll never be religious about this, it'd be great
to provide as many unit tests as we can, while we develop code.

We collect the unit tests in the :repourl:`tests` folder;
:repourl:`tests/test_codestyle.py` is the simplest possible unit
test, while :repourl:`tests/test_spline.py` is an actual working
example. The file names for all the unit-testing python modules should start
with `test_`, because that is the pattern that the test
discovery will look for.

To run the full suite:

.. code-block:: bash

    make test

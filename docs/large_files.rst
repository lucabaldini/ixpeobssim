.. _large_files:

Large data files
================

When setting up simulation models or analyzing data it is sometimes the case that
large data files (say larger than 10 MB) are necessary, XPSEC table models
and Chandra event lists to be fed into ixpeobssim being typical examples.

.. warning::

  Please refrain from pushing large data files to the ixpeobssim repository---a
  code repository is for code, not for data! Even setting aside the fact that
  data files are essentially impossible to diff (i.e., if and when you push a
  new version of a file you are effectively pushing the entire thing over and
  over again, not just the difference with respect to the previous version), if
  you push 1 GB of FITS files to the repository, you are effectively forcing
  everybody else to download 1 GB of data at the next pull, and you are
  increasing the dimensions of the bare repository in a non-reversible way.

  As a rule of thumb, try and avoid pushing files larger than a few MB, unless
  there is a compelling reason for doing so. (Really, once you have pushed a
  file there's essentially no way to get rid of it in git.)

If you have large auxiliary files that are potentially useful for other collaborators there is a
`separate repository <https://bitbucket.org/ixpesw/ixpeobssim_auxfiles/src/master/>`_
that we (ab)use specifically for this purpose, so that only the subset of users
that care about a specific set of files needs to worry about. You will notice that we
don't actually *push* files on that repository---we attach files in the associated
`download page <https://bitbucket.org/ixpesw/ixpeobssim_auxfiles/downloads/>`_, instead.
You should refer to the README file on the
`ixpeobssim_auxfile <https://bitbucket.org/ixpesw/ixpeobssim_auxfiles/src/master/>`_
repository for minimal instructions about the usage of that facility.

On the ``ixpeobssim`` side, as of version 14.0.0, we have a mechanism in place to
direct the user to this repository when a functionality that needs one or more
auxiliary file(s) is used for the first time. Typically all you have to do is
download the necessary file(s) into your ``$IXPEOBSSIM_AUXFILES`` directory
(defaulting to ``~/ixpeobssimauxfiles``), and there are facilities to provide
instructions to the user to do so.
To this end, the :class:`ixpeobssim.srcmodel.magnetar.xMagnetarModelsT2020` interface to
the magnetar table models described in
`Taverna et al. 2020 <https://ui.adsabs.harvard.edu/link_gateway/2020MNRAS.492.5057T/EPRINT_PDF>`_
is the original reason for setting up the large-file infrastructure, and a
good source of inspiration for the associated usage.

In a nutshell, ``ixpeobssim`` provides three convenience functions to handle
auxiliary (large) files, and you can protect the sensitive part of your
configuration or analysis code by a simple call along the lines of

.. code-block:: python

  from ixpeobsim.utils.logging_ import abort
  from ixpeobssim import auxfile_missing

  if auxfiles_missing(*file_list):
      abort()

This should interrupt the execution of the program and print on the terminal
enough information for the user to be able to download the right file(s) in the
right location.

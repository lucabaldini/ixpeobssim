.. _pipeline:

Analysis pipelines
==================

One of the ixpeobssim design goals since the very beginning of the code
development was to allow the user to develop simulation and analysis
pipelines with the minimum possible effort.

The basic technology for developing pipelines is the
:mod:`ixpeobssim.core.pipeline` module.


Application wrappers
--------------------

Each ixpeobssim application, see the :ref:`reference` page, is wrapped in the
:mod:`ixpeobssim.core.pipeline` module so that it can be called from within
a generic Python script with the exact same arguments that one would
pass from command line.

Therefore the shell command 

.. code-block:: bash

    xpobssim --configfile config/toy_point_source.py --duration 10000

should in principle be `exactly` equivalent to the Python snippet

.. code-block:: python

    import ixpeobssim.core.pipeline as pipeline

    file_list = pipeline.xpobssim(configfile='config/toy_point_source.py', duration=10000)

By design all the command-line switches that can be passed when any of the
ixpeobssim applications is invoked in the shell can be passed as keyword
arguments to the corresponding Python wrapper in
:mod:`ixpeobssim.core.pipeline`. Easy, isn't it?
 

Chaining application calls
--------------------------

ixpeobssim Python wrappers typically return the list of all the files that the
function call has created (or, more precisely, a list of strings, each one
representing the path to a specific output file). This makes it easy to chain
applications one after the other, which is a typical use case:

.. code-block:: python

    import ixpeobssim.core.pipeline as pipeline

    evt_file_list = pipeline.xpobssim(configfile='config/toy_disk.py', duration=10000)
    cmap_file_list = pipeline.xpbin(*evt_file_list, algorithm='CMAP')
    count_map = xBinnedMap.from_file_list(cmap_file_list)
    count_map.plot()

While this works as advertised, when dealing with more complex examples
it is often the case that one might want to run specific pieces of the pipeline
independently from the others. In this case passing around file lists with the
mechanism shown in the example above doesn't really make sense, as each step
relies on the fact that all the previous ones did run through.

An alternative possibility that the ixpeobssim pipeline framework offers is
based on the fact that file paths are typically constructed by adding a suffix
to a base name, which in turn coincides with the name of the source model
being simluated and anlyzed.

.. code-block:: python

    import ixpeobssim.core.pipeline as pipeline

    pipeline.setup(model='toy_disk')

    pipeline.xpobssim(configfile='config/toy_disk.py', duration=10000)

    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP')

    file_list = pipeline.file_list('cmap')
    xBinnedMap.from_file_list(file_list)
    count_map.plot()

The reader is referred to the
:repourl:`ixpeobssim/examples/toy_periodic_source.py` example for a somewhat
advanced illustration of the file-list mechanism implemented in the
ixpeobssim pipeline.


Pipeline rc parameters
----------------------

The ixpeobssim pipeline is implemented as a series of stand-alone methods,
and its state is controlled by a look-up dictionary of global parameters that
can be effectively used to exchange informations.

Run-commands parameters are set with the ``setup()`` method and retrieved via
the ``param()`` method

.. code-block:: python

    import ixpeobssim.core.pipeline as pipeline

    pipeline.setup(model=toy_disk)
    print(pipeline.param('model'))

The ``model`` rc param plays a peculiar role, in that under normal conditions
it can be used as a helper to resolve and create file paths. Once the ``model``
parameter is specied one can run, e.g., ixpeobssim without specifying the path
to the configuration file, assuming you want to use the one in the default
folder. The example above can therefore be recasted as

.. code-block:: python

    import ixpeobssim.core.pipeline as pipeline

    pipeline.setup(model='toy_disk')

    pipeline.xpobssim(duration=10000)

    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP')

    file_list = pipeline.file_list('cmap')
    xBinnedMap.from_file_list(file_list)
    count_map.plot()


Pipelines in action
-------------------

The :mod:`ixpeobssim.core.pipeline` module provides a bootstrap function
that should be the preferred way, in practice, to create a simulation and
analysis pipeline. A minimal pipeline example will typically look like


.. code-block:: python

    import ixpeobssim.core.pipeline as pipeline

    def run():
       """Do something useful.
       """
       # Put implementation here
       pass

    if __name__ == '__main__':
       pipeline.bootstrap_pipeline('toy_model')


In this case the boostrap function set the model name for the pipeline,
create a custom option parser that allows to control from command line the
relevant options, and parse the command-line options.

You can run any of the pipelines in the :repourl:`ixpeobssim/examples` folder
with the ``--help`` option to see what the boostrap function makes available,
but in a nutshell this will allow you to

* execute a specific method in your pipeline definition;
* save the plots to file;
* run in batch.

Note that the bootstrap function includes call to the routines showing the
plots and saving them---i.e., you should not call the plt.save() or plt.show()
methods explicitely.

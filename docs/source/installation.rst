Installation & Dependencies
===========================

Both libraries are built using CMake_.

.. _CMake: https://cmake.org/

Dependencies - C++ Libraries
----------------------------

The C++ parts of the code require the following to be installed

- A C++11 compliant compiler
- The GNU Scientific Library (GSL_)
- The ``program_options`` component from the BOOST_ library
- The ROOT_ framework
- The spdlog_ logging functionality

.. _GSL: https://www.gnu.org/software/gsl/
.. _BOOST: http://www.boost.org/doc/libs/1_66_0/doc/html/program_options.html
.. _ROOT: https://root.cern.ch/
.. _spdlog: https://github.com/gabime/spdlog

Additional notes
++++++++++++++++

In Ubuntu the ``program_options`` comes packaged separately and can be installed using

.. code-block:: bash

   sudo apt-get install libboost-program-options-dev

Installation was tested with ROOT version 6, though no problems should occur with previous versions

Dependencies - Python visualization
-----------------------------------

In order to use the Python GUI, one requires

- Python v2.7
- The NumPy_ library
- The Qt4_ framework
- PySide_ python-Qt bindings
- PyQtGraph_
- shell_

.. _NumPy: http://www.numpy.org/
.. _Qt4: http://doc.qt.io/archives/qt-4.8/
.. _PySide: http://wiki.qt.io/PySide
.. _PyQtGraph: http://www.pyqtgraph.org/
.. _shell: https://pypi.python.org/pypi/shell/1.0.1

all of these can be installed using ``pip``

.. code-block:: bash

   sudo pip install numpy pyside pyqtgraph shell

Installation
------------

The libraries do not allow an in-source build. Please run cmake from a separate build folder. Installing should be as simple as running

.. code-block:: bash

   mkdir build && cd build
   cmake ../

All available options can be revised using the graphical ``ccmake`` program. In case the ``spdlog.h`` file was not found, the directory in which it resides can be added to the ``PATH`` environment variable. Running cmake again should resolve the issue. After this one only has to build all targets using

.. code-block:: bash

   make

Installation is optional, and can be run using

.. code-block:: bash

   make install

This may require sudo privileges.

Compiled libraries and executables can be found in the corresponding directories.

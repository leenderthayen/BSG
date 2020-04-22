Installation & Basic execution
******************************

Installation & Dependencies
===========================

Both libraries are built using CMake_.

Databases used in the graphial user interface can be downloaded automatically through the CMake script. See the options.

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

- In Ubuntu the ``program_options`` comes packaged separately and can be installed using

.. code-block:: bash

   sudo apt-get install libboost-program-options-dev

- Installation was tested with ROOT version 6, though no problems should occur with previous versions

- Installing spdlog through your package manager may install an outdated version. Please use the source code on github.

On macOS, you can do that simply by using brew_:

.. _brew: https://brew.sh/

.. code-block:: bash

   brew install root spdlog boost


Dependencies - Python visualization
-----------------------------------

In order to use the Python GUI, one requires

- Python v2.7
- The NumPy_ library
- The Qt4_ framework
- PySide_ python-Qt bindings
- PyQtGraph_
- shell_
- The QDarkStyle_ style sheet (optional)

.. _NumPy: http://www.numpy.org/
.. _Qt4: http://doc.qt.io/archives/qt-4.8/
.. _PySide: http://wiki.qt.io/PySide
.. _PyQtGraph: http://www.pyqtgraph.org/
.. _shell: https://pypi.python.org/pypi/shell/1.0.1
.. _QDarkStyle: https://github.com/ColinDuquesnoy/QDarkStyleSheet

all of these can be installed using ``pip`` using the requirements file provided. All of this is performed automatically when using CMake.

.. code-block:: bash

   sudo pip install -r requirements.txt

Note that support for Qt5 is possible, using PySide 2 to convert the MainWindow.ui file into Qt5-friendly code, and installing the github version of pyqtgraph.

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

Compiled libraries and binaries can be found in the corresponding directories.

For smooth use of the GUI, the following paths can be defined in your ``.bashrc``:

- ``BSGPATH``: absolute path to the compiled ``bsg_exec`` executable.
- ``BSGEXCHANGEPATH``: absolute path to the ``ExchangeData.dat`` file.
- ``ENSDFDIR``: absolute path to the full ENSDF_ library, downloaded using cmake if the option is enabled
- ``FRDMPATH``: absolute path to the FRDM_ library, downloaded using cmake if the option is enabled
- ``CHARGERADIIPATH``: absolute path to the ChargeRadii_ library, downloaded using cmake if the option is enabled

.. _ENSDF: https://www.nndc.bnl.gov/ensarchivals/
.. _FRDM: https://www.sciencedirect.com/science/article/pii/S0092640X1600005X
.. _ChargeRadii: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.94.064315

Execution
=========

Execution of the program is performed as with any other linux program. Make sure to either specify the location of the config.txt and ExchangeData.dat files or the location, or have them (using a soft link) in the current folder.

Using the 63Ni beta decay as an example, executation could as simple as

.. code-block:: bash

   ./bsg_exec -i 63Ni.ini -o 63Ni
                
upon which 4 files will be created detailing the calculation. The file ending in .txt contains a general overview. Note, *you* have to create the file 63Ni.ini, we'll see how you do that later in the next sections.

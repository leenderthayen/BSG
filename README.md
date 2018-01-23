BSG -- Beta Spectrum Generator
==============================
![alt text](https://img.shields.io/badge/License-MIT-blue.svg 'License')
![alt text](https://img.shields.io/badge/Python-2.7-blue.svg 'Python version')
![alt text](https://img.shields.io/badge/Linux-Supported-brightgreen.svg 'Supported OS')
[![Documentation Status](https://readthedocs.org/projects/bsg/badge/?version=latest)](http://bsg.readthedocs.io/en/latest/?badge=latest)

Purpose
-------
Contributors:
* Leendert Hayen (leendert.hayen@kuleuven.be)

This package calculates nuclear matrix elements relevant to allowed beta decay, and its corresponding spectrum according to the formalism by Hayen *et al.* ([arXiv](https://arXiv.org/abs/1709.07530)). It contains a graphical user interface written in Python which allows connection to, e.g., the ENSDF database.

Dependencies
------------
The C++ component of this package makes use of the following libraries:
* [GSL](https://www.gnu.org/software/gsl/)
* [boost::program_options](http://www.boost.org/doc/libs/1_66_0/doc/html/program_options.html)
* [ROOT](https://root.cern.ch/)
* [spdlog](https://github.com/gabime/spdlog)

Installation requires C++11.

The Python GUI makes use of the following libraries:
* [NumPy](http://www.numpy.org/)
* [Qt4](http://doc.qt.io/archives/qt-4.8/)
* [PySide](http://wiki.qt.io/PySide)
* [PyQtGraph](http://www.pyqtgraph.org/)
* [shell](https://pypi.python.org/pypi/shell/1.0.1)

Installation
------------
A CMake script is provided to automate the installation. Detailed documentation can be found [here](http://bsg.readthedocs.io).

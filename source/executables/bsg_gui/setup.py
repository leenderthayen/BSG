from distutils.core import setup

setup(
    name='BSG_GUI',
    version='0.1.0',
    author='L. Hayen',
    author_email='leendert.hayen@gmail.com',
    packages=['',],
    scripts=['',],
    url='',
    license='LICENSE.txt',
    description='Graphical user interface for the Beta Spectrum Generator library.',
    long_description=open('README.txt').read(),
    requires=[
        "shell (>= 1.0.1)",
        "QDarkStyle (>= 2.6.4)",
	"PySide (>= 1.2.4)",
	"configparser (>= 3.5.0)",
	"numpy"
    ],
)

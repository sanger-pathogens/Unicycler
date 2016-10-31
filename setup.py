#!/usr/bin/env python3
"""
Run 'python3 setup.py install' to install Unicycler.
"""

# Make sure this is being run with Python 3.4 or later.
import sys
if sys.version_info.major != 3 or sys.version_info.minor < 4:
    sys.exit('Error: you must execute setup.py using Python 3.4 or later')

import os
import shutil
import shlex
from distutils.core import Command
import subprocess
import multiprocessing
import fnmatch
import importlib.util

# Install setuptools if not already present.
if not importlib.util.find_spec("setuptools"):
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup
from setuptools.command.install import install

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

# Get the program version from another file.
exec(open('unicycler/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCRIPTION = readme.read().decode()


class UnicycleInstall(install):
    """
    The install process copies the C++ shared library to the install location.
    """
    user_options = install.user_options + [
        ('makeargs=', None, 'Arguments to be given to make')
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.makeargs = None

    def run(self):
        # Make sure we have permission to write the files.
        if os.path.isdir(self.install_lib) and not os.access(self.install_lib, os.W_OK):
            sys.exit('Error: no write permission for ' + self.install_lib + '  ' +
                     'Perhaps you need to use sudo?')
        if os.path.isdir(self.install_scripts) and not os.access(self.install_scripts, os.W_OK):
            sys.exit('Error: no write permission for ' + self.install_scripts + '  ' +
                     'Perhaps you need to use sudo?')

        # Clean up any previous Unicycler compilation.
        clean_cmd = ['make', 'distclean']
        self.execute(lambda: subprocess.call(clean_cmd), [],
                     'Cleaning previous compilation: ' + ' '.join(clean_cmd))

        # Build Unicycler's C++ code.
        make_cmd = ['make']
        try:
            make_cmd += ['-j', str(min(8, multiprocessing.cpu_count()))]
        except NotImplementedError:
            pass
        if self.makeargs:
            make_cmd += shlex.split(self.makeargs)
        self.execute(lambda: subprocess.call(make_cmd), [],
                     'Compiling Unicycler: ' + ' '.join(make_cmd))
        cpp_code = os.path.join('unicycler', 'cpp_functions.so')
        if not os.path.isfile(cpp_code):
            sys.exit("Error: compilation of Unicycler's C++ component failed")

        install.run(self)

        # Copy non-Python stuff to the installation directory.
        shutil.copyfile(cpp_code, os.path.join(self.install_lib, 'unicycler', 'cpp_functions.so'))
        gene_data_source_dir = os.path.join('unicycler', 'gene_data')
        gene_data_dest_dir = os.path.join(self.install_lib, 'unicycler', 'gene_data')
        if not os.path.exists(gene_data_dest_dir):
            os.makedirs(gene_data_dest_dir)
        shutil.copyfile(os.path.join(gene_data_source_dir, 'start_genes.fasta'),
                        os.path.join(gene_data_dest_dir, 'start_genes.fasta'))
        shutil.copyfile(os.path.join(gene_data_source_dir, 'lambda_phage.fasta'),
                        os.path.join(gene_data_dest_dir, 'lambda_phage.fasta'))

        # Display a success message!
        try:
            from unicycler.misc import get_ascii_art
            print(get_ascii_art())
        except ImportError:
            print()
        print('Unicycler is installed!')
        print()
        print('Example commands:')
        print('  unicycler --help')
        print('  unicycler --help_all')
        print('  unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz '
              '-l long_reads.fastq.gz -o path/to/output_dir')
        print('  unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz --no_long '
              '-o path/to/output_dir')
        print()


class UnicycleClean(Command):
    """
    Custom clean command that really cleans up everything, except for:
      - the compiled *.so file needed when running the programs
      - setuptools-*.egg file needed when running this script
    """
    user_options = []

    def initialize_options(self):
        self.cwd = None

    def finalize_options(self):
        self.cwd = os.getcwd()

    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd

        delete_directories = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for dir_name in fnmatch.filter(dir_names, '*.egg-info'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, 'build'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, '__pycache__'):
                delete_directories.append(os.path.join(root, dir_name))
        for delete_directory in delete_directories:
            print('Deleting directory:', delete_directory)
            shutil.rmtree(delete_directory)

        delete_files = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for filename in fnmatch.filter(filenames, 'setuptools*.zip'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.o'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.pyc'):
                delete_files.append(os.path.join(root, filename))
        for delete_file in delete_files:
            print('Deleting file:', delete_file)
            os.remove(delete_file)


setup(name='unicycler',
      version=__version__,
      description='bacterial genome assembler for hybrid read sets',
      long_description=LONG_DESCRIPTION,
      url='http://github.com/rrwick/unicycler',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPL',
      packages=['unicycler'],
      entry_points={"console_scripts": ['unicycler = unicycler.unicycler:main',
                                        'unicycler_align = unicycler.unicycler_align:main',
                                        'unicycler_check = unicycler.unicycler_check:main',
                                        'unicycler_polish = unicycler.unicycler_polish:main']},
      zip_safe=False,
      cmdclass={'install': UnicycleInstall,
                'clean': UnicycleClean}
      )

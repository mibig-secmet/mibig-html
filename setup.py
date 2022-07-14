"""Setuptools magic to install MIBiG HTML generator."""
import glob
import os
import subprocess
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


def read(fname):
    """Read a file from the current directory."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


long_description = read('README.md')

install_requires = [
    'antismash @ git+https://github.com/antismash/antismash@6-1-1',
    'eutils',
    'python-mibig',
    'numpy',
    'mibig-taxa',
    'requests',
]

tests_require = [
    'pytest >= 3.4.0, < 5',  # pytest 5 breaks compatibility with coverage
    'coverage',
    'pylint == 2.10.2',
    'mypy == 0.910',  # for consistent type checking
    'types-requests',
]


def read_version() -> str:
    """Read the version fromt he appropriate place in the library."""
    with open(os.path.join('mibig_html', 'main.py'), 'r') as handle:
        for line in handle:
            if line.startswith('__version__'):
                return line.split('=')[-1].strip().strip('"')
    return "unknown"


def find_data_files():
    """Setuptools package_data globbing is stupid, so make this work ourselves."""
    data_files = []
    for pathname in glob.glob("mibig_html/**/*", recursive=True):
        if pathname.endswith('.pyc'):
            continue
        if pathname.endswith('.py'):
            continue
        if '__pycache__' in pathname:
            continue
        pathname = glob.escape(pathname)
        pathname = pathname[10:]
        data_files.append(pathname)
    if "HARDCODE_MIBIG_GIT_VERSION" in os.environ:
        version_file = os.path.join('mibig_html', 'git_hash')
        with open(version_file, 'wt') as handle:
            try:
                git_version = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                                      universal_newlines=True).strip()
                changes = subprocess.check_output(['git', 'status', '--porcelain'],
                                                  universal_newlines=True).splitlines()
                if len(changes) != 0:
                    git_version += "(changed)"
                handle.write(git_version)
            except (OSError, subprocess.CalledProcessError):
                pass
        data_files.append(version_file)
    return data_files


class PyTest(TestCommand):
    """Allow running tests via python setup.py test."""

    def finalize_options(self):
        """Test command magic."""
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        """Run tests."""
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name="mibig_html",
    python_requires='>=3.7',
    version=read_version(),
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_data={
        'mibig_html': find_data_files(),
    },
    author='antiSMASH development team',
    author_email='antismash@secondarymetabolites.org',
    description='The antibiotics and Secondary Metabolites Analysis Shell.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_requires,
    tests_require=tests_require,
    entry_points={
        'console_scripts': [
            'mibig_html=mibig_html.__main__:entrypoint',
        ],
    },
    cmdclass={'test': PyTest},
    url='https://github.com/mibig-secmet/antismash-mibig',
    license='GNU Affero General Public License v3 or later (AGPLv3+)',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Operating System :: OS Independent',
    ],
    extras_require={
        'testing': tests_require,
    },
)

import os
from setuptools import setup

# read in README
this_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_dir, 'README.md'), 'rb') as f:
    long_description = f.read().decode().strip()

# requirements
install_requires = [
    "spiceypy",
    "numpy",
    "pandas",
    "matplotlib"
]
extras_require = {
    'dev': [
        'autopep8',
        'flake8',
        'pytest >= 5.0, < 5.4',
        'pytest-console-scripts',
        'pytest-cov',
        'pytest-runner',
        'twine',
    ],
    'docs': [
        'sphinx',
        'sphinx_rtd_theme',
        'sphinxcontrib-programoutput',
    ],
}

setup(
    name = 'oif',
    description = 'A package to efficiently tell whether a set of objects are in field',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'http://github.com/AsteroidSurveySimulator/objectsInField',
    author = 'Shantanu Naidu and the OIF Team',
    author_email = 'shantanu.p.naidu@jpl.nasa.gov',
    license = 'BSD 3-Clause',

    packages = ['oif'],

    entry_points = {
        'console_scripts': [
            'oif = oif.__main__:main',
        ],
    },
    package_data = {
        "oif": ["data/*"]	# this is where various standard files and kernels live
    },
    exclude_package_data = {
        "": ["*~"]		# don't include backups
    },

    python_requires = '>=3.6.*',
    install_requires = install_requires,
    extras_require = extras_require,
    setup_requires = ['setuptools_scm'],
    zip_safe=False,
    use_scm_version = {
        'write_to': 'oif/_version.py'
    },

    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'License :: OSI Approved :: BSD License',
    ],
)

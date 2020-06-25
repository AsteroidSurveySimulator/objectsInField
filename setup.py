import os
from setuptools import setup

##
# Infer version from git tag information
##

def git_version():
    # Uses the following logic:
    #
    # * If the current commit has an annotated tags, the version is simply the tag with
    #   the leading 'v' removed.
    #
    # * If the current commit is past an annotated tag, the version is constructed at:
    #
    #    '{tag+1}.dev{commitcount}+{gitsha}'
    #
    #  where {commitcount} is the number of commits after the tag (obtained with `git describe`)
    #  and {tag+1} is the tag version incremented by one (i.e., if tag=0.0.1, it's 0.0.2). This
    #  is needed because .dev tags compare as less than non-.dev tags, per PEP-440.
    #
    # * If there are no annotated tags in the past, the version is:
    #
    #    '0.0.0.dev{commitcount}+{gitsha}'
    #
    # Inspired by https://github.com/pyfidelity/setuptools-git-version
    # Creates PEP-440 compliant versions

    from subprocess import check_output

    command = 'git describe --tags --long --dirty --always'
    version = check_output(command.split()).decode('utf-8').strip()

    parts = version.split('-')
    if len(parts) in (3, 4):
        dirty = len(parts) == 4
        tag, count, sha = parts[:3]
        if not tag.startswith('v'):
            raise Exception(
                "Annotated tags on the repository must begin with the letter 'v'. Please fix this then try building agains.")
        tag = tag[1:]
        if count == '0' and not dirty:
            return tag
    elif len(parts) in (1, 2):
        tag = "0.0.0"
        dirty = len(parts) == 2
        sha = parts[0]
        # Number of commits since the beginning of the current branch
        count = check_output(
            "git rev-list --count HEAD".split()).decode('utf-8').strip()

    # this will be a .dev version; assume the version is of the form
    # major[.minor[.patchlevel]], and increment the patchlevel by 1.
    # because minor and/or patchlevel are optional
    (major, minor, patchlevel) = (tag + ".0.0").split('.')[:3]
    patchlevel = str(int(patchlevel) + 1)
    tag = '.'.join((major, minor, patchlevel))

    # if the working directory is dirty, append a '.dYYYYMMDD' tag
    if dirty:
        import time
        dirtytag = time.strftime('.d%Y%m%d')
    else:
        dirtytag = ''

    fmt = '{tag}.dev{commitcount}+g{gitsha}{dirtytag}'
    return fmt.format(tag=tag, commitcount=count, gitsha=sha.lstrip('g'), dirtytag=dirtytag)

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
    version=git_version(),
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
    zip_safe=False,

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

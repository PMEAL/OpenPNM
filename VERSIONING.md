## Release Management and Versioning

`OpenPNM` uses [Semantic Versioning](http://semver.org) (i.e. X.Y.Z) to label releases.  All major and minor versions (X.Y.z) are available on [PyPI](https://pypi.python.org/pypi), but bugfix releases (x.y.Z) are not generally pushed unless the bug is important.  Bugfix releases are available via download of the source code from Github.

`OpenPNM` uses the [Github Flow](https://guides.github.com/introduction/flow/) system of Git branching, except instead of merging PRs into *master*, they are merged into a branch called *dev*. Any code added to *dev* is done via Pull Requests (PRs).  When new PRs are merged into the *dev* branch, they are *not* given a new version number. Once enough new features have been added, the *dev* branch is merged into the *master* branch, and the minor release number (x.Y.z) will be incremented. An exception to this rule are bugfixes which may be found on *master*.  In these cases a PR can be merged into *master* and the version number will be incremented (x.y.Z) to indicate the fix.

`OpenPNM` depends on several other packages widely known as the [Scipy Stack](https://www.scipy.org/stackspec.html).  It is our policy to always support the latest version of all these packages and their dependencies.

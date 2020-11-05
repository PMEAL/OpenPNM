## Release Management and Versioning

`OpenPNM` uses [Semantic Versioning](http://semver.org) (i.e. X.Y.Z) to label releases.  As of verion 2.3.2, all versions of OpenPNM are available on [PyPI](https://pypi.python.org/pypi).  Prior to this, only major and minor version were pushed.

All development occurs on `dev` via feature branches and the pull request functionality of Github. A new release is defined each time the `dev` branch is merged into the `release` branch. Several automations are setup so that upon each release, the code is automatically deployed to PyPi and Conda, and a release announcement is created on Github containing a summary of all the changes.  This `dev` and `release` workflow replaces the previous approach based on gitflow.

`OpenPNM` depends on several other packages widely known as the [Scipy Stack](https://www.scipy.org/stackspec.html).  It is our policy to always support the latest version of all these packages and their dependencies.

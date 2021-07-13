# Running development version of OpenPNM

The OpenPNM development process consists of a ``release`` branch, which gets deployed to *conda* and *pip* and is available via ``conda install -c conda-forge openpnm`` and ``pip install openpnm``, respectively; while development of new features occurs on the ``dev`` branch.  We try to regularly merge ``dev`` into ``release``, but sometimes several weeks or even months can elapse between releases. If you really need a feature that is only available on the ``dev`` branch, OR you want to fiddle with the source code yourself (hopefully with the intention of making a contribution!), you need to install OpenPNM in a different way. This will be explained in this example.

## Install a git client

To begin, you will need the ability to interact with a "git repository" ([Git is a version control system](https://en.wikipedia.org/wiki/Git) that is essential for software development).  This *can* be done by installing the [official git package](https://git-scm.com/), but this is a command-line interface and requires a fair bit of expertise.  

Luckily there are numerous high-quality graphical interfaces for git, including:

* [sourcetree](https://www.sourcetreeapp.com/)
* [fork](https://git-fork.com/)
* [gitkraken](https://www.gitkraken.com/)

to name just a few of our favorites.  Each has some limitations.  For instance GitKraken works on all the major OS's, but is not free for private repos, while Fork works only on Mac and Windows. Choose the one that's right for you and install it.

## Clone the OpenPNM repo from Github

You now need to download the code github onto your computer, which can be done through your git client.  Each git client is sligthly different, but the following steps are fairly general:

* There should be an option to *clone* a repo, which asks for the location of the remote files (https://github.com/PMEAL/OpenPNM.git), and the location where you'd like the local copy to be stored (e.g. C:\users\yourname\code\OpenPNM).  
* Within your git client you can select which branch and version of OpenPNM you would like to use (probably by double-clicking the branch name), this is called a *checkout*. 
* By default most git clients will only download the default branch (which is `dev` in OpenPNM) so if you need to *checkout* a different branch, you have to browse under the "origin" which refers to the code on Github. 
* Checking out a branch from the "origin" will make that code available on your local machine.

## Install the local version

Having the code locally is not enough to run it, since you must tell python about this situation, using the following steps:

* If you have installed OpenPNM from either pip or conda, you should uninstall it first (``pip uninstall openpnm`` or ``conda remvoe openpnm``).  
* Next you navigate to the folder where you saved the OpenPNM git repo (e.g. ``cd C:\users\username\code\OpenPNM``) and install using pip with the following:  ``pip install -e .``.  
  * The ``-e`` tells pip that this code is 'editable', meaning it should reload the current version of the files on each import.  This is necessary so changes made the files (like switching branches) are reflected in the package.  
  * The ``.`` tell pip to find the ``setup.py`` file in the "current directory" and run it.  
  * Alternatively you could have told pip where to look, such as ``pip install -e "C:\users\username\code\OpenPNM"``. 

> Note: As of version 2.5 OpenPNM relies on "conda" to install packages since we use some package that cannot be distributed via pip (i.e. pypardiso).  So you may need to run ``conda install -c conda-forge pypardiso``...though you probably already have this if you install openpnm from conda previously.

## Conributing to OpenPNM

The above is sufficient to get you using the ``dev`` version, however, if you are wishing to contribute back to OpenPNM (and we hope you do) you must do a few more steps.  
* Firstly, you should make your changes on a 'branch', not on ``dev``.  You can create branches easily in your git client, then commit your changes to the branch.  
* BUT, you cannot 'push' your code changes to a repository on Github without permission of the owners.  Instead you need to [*fork*](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) our repo into your own Github account. You can then push your new branch to your own version of the repo (since you own it). 
* You must then ask us to "pull" your changes from your repo onto our repo using the ["pull request"](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/about-pull-requests) functionality of Github.  During this process we can discuss your proposed changes, suggest edits, and review the code before (hopefully) accepting your contribution and merging it into the official OpenPNM package.  
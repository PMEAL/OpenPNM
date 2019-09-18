```
THIS REPO IS NO LONGER ACTIVELY MAINTAINED. WE HAVE MOVED ALL THE EXAMPLES ONTO THE [MAIN REPO](https://github.com/PMEAL/OpenPNM/tree/master/examples).  THE REASON FOR DOING THIS IS SO THAT THE EXAMPLES WILL ALL BE TESTED WITH EACH NEW VERSION OF THE CODE TO ENSURE THEY WORK. 
```

[![Travis](https://travis-ci.org/PMEAL/OpenPNM.svg?branch=master)](https://travis-ci.org/PMEAL/OpenPNM)

# OpenPNM-Examples

A collection of scripts illustrating the use of OpenPNM.  The files in these folders are written in MarkDown, which allows for an useful blend of formatted text and highlighted code.  You can browse the available files by clicking through the web interface of the repository seen above.
> **Note:** All the examples in this repository are subject to automated testing, so if the build badge above says passing then all of these examples should work.  All code within these files is in code-blocks and preceded by the python prompt symbol (>>>).  This is annoying for cutting-and-pasting, but is necessary if we want to run automated tests on the code.

## Contents

The examples are organized in the directory structure of the repository.  The folders and files are given somewhat descriptive names to indicate what sort of example they contain.

## How to Contribute an Example to this Collection

We strongly encourage you to create an example and post it in this collection.  The procedure is straightforward and can be explained in just a few lines:

1.  Fork this repository to your own Github account using the 'Fork' button on the top of this page.  This will create an independent version of the repo under your own Github account, and will let you make edits and additions to the repository without needing write permission from us.

2. Create your example using the instructions below, and push your changes to your forked version of the repository.  You can edit and push as many times as needed to check how the document looks when Github renders it.

3. Make a pull request using the Github website by going the 'Pull Requests' tab on your forked repository, and then selecting the 'New pull request' button.  This will trigger the automated testing to ensure your example works.  You can help ensure this will pass by running the tests locally by typing `python run_tests` at your command prompt (assuming you have [Pytest installed](https://pytest.org/latest/getting-started.html)).

## How to Write an Example

This section introduces a few rules and tips for creating nicely formatted documents that would be suitable for inclusion here.

1. All files must be written in MarkDown.  This is the text formatting language used by Github to create nice looking issues and posts within a given repository.  MarkDown files have an *md* file extension.  If you need some help on MarkDown formatting and notation, checkout the [MarkDown Cheatsheet].  The style is not very strict, but efforts should be made to be consistent with the existing files.

2. All code must be not only enclosed with code blocks denoted by triple ticks (\`\`\`)  (see the [MarkDown Cheatsheet]), but each line must be preceded by >>> like the Python prompt.  The Pytest package which is used for automated testing of examples only recognizes lines of code preceded by the a >>>.  Moreover, each code block must end with one blank line before the closing ticks (\`\`\`), or else Pytest reports an error.  The flip side of this coin is that if you feel that some could should NOT be tested, this can be accomplished by excluding the >>>.

3. Any images to be inserted into an examples should be hosted on an external size such as [imgur](http://imgur.com).  These can then be linked to using the appropriate notations as given in the [MarkDown Cheatsheet].

4. Any files required by an example (such as data files containing network info) should be stored in the 'fixtures' directory at the top level of this repository.

5.  The file name should be descriptive of the contents since the user only sees these when navigating the repository. Similarly, the file should be placed in an appropriate folder, or a folder should be created if necessary.

[MarkDown CheatSheet]: https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet

# Developer's guide

## How to submit an [issue](https://github.com/PMEAL/OpenPNM/issues)

In case any of the following applies, go ahead and create an [issue](https://github.com/PMEAL/OpenPNM/issues) *.

- You have a general question such as "how to do X in `OpenPNM`" **
- You have found a bug or even a typo in our documentation or in the dosctrings
- You have a nice or essential feature in mind that's missing in `OpenPNM`

\* When submitting a new issue, we ask that you use one of the three available templates. The templates have sections that give sufficient information and context to the developers to address your issue as quickly as possible.

\** We also insist that you first go through our [examples](https://github.com/PMEAL/OpenPNM/tree/dev/examples). Our examples are fairly well organized, and there is a table of contents that you should be able to skim through fairly quickly.

## How to submit a [pull request](https://github.com/PMEAL/OpenPNM/pulls)

When proposing a change of any sort (bug fix, new feature, typo in documentation, etc.), it is recommended that you first open up an [issue](https://github.com/PMEAL/OpenPNM/issues), explaining the bug, new feature, etc, and discussion with the developers can occur to help focus the issue.

Once the issue and solution are well defined, you will hopefully proceed with creating the pull request to actually implement the code changes.  The steps for creating a pull request are outlined in [this tutorial](https://github.com/PMEAL/OpenPNM/blob/dev/examples/tutorials/using_and_contributing_to_the_dev_version.md).

GitHub has a way of associating issues with pull requests (see [here](https://docs.github.com/en/free-pro-team@latest/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue#linking-a-pull-request-to-an-issue-using-a-keyword) for more details). You can use the following keywords to do this:

- `close` , `closes` , `closed`
- `fix` , `fixes` , `fixed`
- `resolve` , `resolves` , `resolved`

For instance, if your pull request is addressing issues 15, 28, and 122, you can include the following message somewhere in the pull request:

```
Fixes #15, closes #28, resolves #122.
```

The important part is that you need to use multiple keywords for multiple issues (you can use the same keyword though, e.g. using `fixes` multiple times).

## Development workflow

We use semantic versioning for versioning `OpenPNM`. In case you're not familiar, it's basically making version numbers sort of human-readable! A version number must be in form of three digits separated by a dot, for example `3.2.5`. You only change the first digit whenever you make a "major" change to your software, for instance changes that make your software backwards-incompatible. The second digit is changed when you make a "minor" change, and finally for bug fixes and other "patch"-type changes, you bump the third digit.

The way we develop `OpenPNM` is as follows:

- We always merge pull requests onto the `dev` branch. So, you can get the latest and greatest changes on the `dev` branch. Steps for cloning the repo using a git client are outlined in [this tutorial](https://github.com/PMEAL/OpenPNM/blob/dev/examples/tutorials/using_and_contributing_to_the_dev_version.md)
- Once in a while (usually three to four times a year), we merge the `dev` branch onto the `release` branch. As the name suggests, the `release` branch is the one that is installed when using `pip` or `conda` (don't use `pip`, use `conda` instead :sunglasses:). Although we don't merge any broken pull request onto the `dev` branch we often tweak things before final release, so generally the `release` branch is more stable.
- When merging the `dev` branch onto the `release` branch, our automated workflow is triggered, generating the changelog and publishing `OpenPNM` on `PyPI` and `conda`. Here's an overview of what's happening under the hood for those interested!

<p align="center">
<img src="https://user-images.githubusercontent.com/14086031/97095215-a9e92b00-162a-11eb-89b7-90c2de650822.png" width=60%>
</p>

## Merging pull requests

For our changelog generator to work well, the only thing you need to remember is to use certain "keywords" when merging onto the `dev` and the `release` branch.

### Merging onto `dev`

When merging your branch onto `dev`, as the merge message, describe what your pull request does concisely and preferably in a single sentence plus one of the following "standard" keywords:

| Change type  | Standard keywords *  | Magical keywords **                                   |
|:-------------|:---------------------|:------------------------------------------------------|
| New feature  | `#new`               | `feature` , `added`                                   |
| Enhancement  | `#enh`               | `revamped` , `improved` , `enhanced` , `optimized`    |
| Maintenance  | `#maint`             | `backend`                                             |
| API change   | `#api`               | `deprecated` , `changed` , `removed` , `modified`     |
| Bug fix      | `#bug`               | `bugfix` , `hotfix` , `fixed`                         |
| Documentation| `#doc`               | `documentation` , `docstring`                         |

\* **Standard keywords**: For consistency, make sure you always use these. Example merge commit message:

```
`topotools.plot_connections` no longer accepts list of pores [#api].
```

\** **Magical keywords**: For ease of use - and also in case you forget to use standard keywords -, feel free to use these in your merge commit message, they will automatically get caught. Example merge commit message:

```
Optimized `topotools.find_neighbor_pores` which is now ~10X faster.
```

**Note 1**: You can use multiple keywords in case your pull request is doing multiple things (e.g. fixes a bug + deprecates a method), although this is discouraged. Please make a separate pull request for every change.

**Note 2**: We're deprecating the magical keywords from `v2.6.0`, so only use the special keywords.

### Merging `dev` onto `release`
Finally, if we're ready to publish a new release to PyPI and `conda`, you should create a pull request, asking to merge the `dev` branch onto the `release` branch. Again, this process is automated so that the version number gets bumped accordingly.  The only thing you need to remember is to use the proper keyword, so that our automated workflow knows how to bump the version number. Please use the following keywords:

- For patches, e.g. bug fixes, etc. use `#patch`
- For minor changes, e.g. a new feature, or a non-breaking API change, use `#minor`
- For major changes, e.g. major backwards-incompatible changes, use `#major`

An example merge commit message:

```
New minor release [#minor]
```



import git
import OpenPNM as op

repo = git.Repo(search_parent_directories=True)
assert op.__version__ == repo.git.describe("--tags").strip('vV')

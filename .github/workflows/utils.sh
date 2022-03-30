# -------------------------------------------------------------------------- #
# Retrieves the most recent tag of the git project in the current working
# directory.
#
# Notes
# -----
# This method is equivalent to calling `get_nth_recent_tag 1`
#
# Example
# -------
# For instance, suppose a project has the following tags: v1.2, v2.0, v3.5
# - get_most_recent_tag
# >>> v3.5
# -------------------------------------------------------------------------- #
function get_most_recent_tag {
    local temp=$(get_nth_recent_tag 1)
    echo "$temp"
}


# -------------------------------------------------------------------------- #
# Retrieves the version number given the location of the version file
#
# Parameters
# ----------
# version_loc: location of the file containing the version number
#
# Notes
# -----
# This method assumes the global variable 'version_loc' holds the relative path
# to the version file.
#
# Example
# -------
# - bump_version minor
# - bump_version patch
# -------------------------------------------------------------------------- #
function get_version {
    local version_loc=$1
    local temp=$(grep -E -o "([0-9]{1,}\.)+[0-9]{1,}(.dev[0-9]{1,})?" $version_loc)
    echo "$temp"
}


# -------------------------------------------------------------------------- #
# Retrieves the nth most recent tag of the git project in the current working
# directory.
#
# Parameters
# ----------
# n: nth recent tag
#
# Example
# -------
# For instance, suppose a project has the following tags: v1.2, v2.0, v3.5
# - get_nth_recent_tag 1
# >>> v3.5
# - get_nth_recent_tag 3
# >>> v1.2
# -------------------------------------------------------------------------- #
function get_nth_recent_tag {
    git fetch --all --tags --force >/dev/null
    local tags=($(git for-each-ref --sort=-creatordate --format '%(refname:strip=2)' refs/tags --count=$1))
    echo "${tags[$(($1-1))]}"
}


# -------------------------------------------------------------------------- #
# Bumps version number based on the release type parameter that could be
# 'patch', 'minor', or 'major'.
#
# Parameters
# ----------
# bump_type: release type, choose from ['patch', 'minor', 'major']
# version_loc: location of the file containing the version number
#
# Notes
# -----
# This method assumes the global variable 'version_loc' holds the relative path
# to the version file.
#
# Example
# -------
# - bump_version minor openpnm/__version__.py
# - bump_version patch openpnm/__version__.py
# -------------------------------------------------------------------------- #
function bump_version {
    local bump_type=$1
    local version_loc=$2
    local version=$(get_version $version_loc)
    bump2version --current-version $version $bump_type $version_loc --verbose
}


# -- Test cases -- #
# git checkout .
# version_loc="openpnm/__version__.py"
# get_most_recent_tag
# get_nth_recent_tag 2
# get_version $version_loc
# bump_version patch $version_loc

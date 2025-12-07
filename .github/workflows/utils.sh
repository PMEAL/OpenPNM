#!/usr/bin/env bash
# Utility functions for version and tag management.

_log_utils() {
    echo "[utils] $*" >&2
}

# -------------------------------------------------------------------------- #
# Retrieves the most recent tag of the git project in the current working directory.
# Usage: get_most_recent_tag
# -------------------------------------------------------------------------- #
get_most_recent_tag() {
    _log_utils "Getting most recent tag..."
    get_nth_recent_tag 1
}

# -------------------------------------------------------------------------- #
# Retrieves the version number from the provided file location.
# Usage: get_version <version_file_path>
# -------------------------------------------------------------------------- #
get_version() {
    local version_loc="$1"
    _log_utils "Reading version from: $version_loc"
    if [[ -z "$version_loc" ]]; then
        echo "Error: version file path not provided" >&2
        return 1
    fi
    if [[ ! -f "$version_loc" ]]; then
        echo "Error: version file not found at $version_loc" >&2
        return 1
    fi
    local version
    version=$(grep -E -o "([0-9]{1,}\.)+[0-9]{1,}(.dev[0-9]{1,})?" "$version_loc" | head -n1)
    _log_utils "Found version: $version"
    echo "$version"
}

# -------------------------------------------------------------------------- #
# Retrieves the nth most recent tag of the git project.
# Usage: get_nth_recent_tag <n>
# -------------------------------------------------------------------------- #
get_nth_recent_tag() {
    local n="$1"
    _log_utils "Fetching tag #$n (1=most recent)"
    if ! [[ "$n" =~ ^[0-9]+$ ]]; then
        echo "Error: Argument must be a positive integer" >&2
        return 1
    fi
    _log_utils "Fetching tags from remote..."
    git fetch --tags --force --quiet
    local tags
    read -ra tags <<< "$(git for-each-ref --sort=-creatordate --format '%(refname:strip=2)' refs/tags --count="$n" | tr '\n' ' ')"
    _log_utils "Found ${#tags[@]} tag(s): ${tags[*]}"
    if (( ${#tags[@]} < n )); then
        echo "Error: Less than $n tags found (found ${#tags[@]})" >&2
        return 1
    fi
    local result="${tags[$((n-1))]}"
    _log_utils "Returning tag: $result"
    echo "$result"
}

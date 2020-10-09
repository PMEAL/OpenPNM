function get_most_recent_tag {
    temp=$(get_nth_recent_tag 1)
    echo "$temp"
}


function get_version {
    version_loc=$1
    temp=$(egrep -o "([0-9]{1,}\.)+[0-9]{1,}" $version_loc/__version__.py)
    echo "$temp"
}


function get_nth_recent_tag {
    git fetch --all --tags
    tags=($(git for-each-ref --sort=-creatordate --format '%(refname:strip=2)' refs/tags --count=$1))
    echo "${tags[$(($1-1))]}"
}

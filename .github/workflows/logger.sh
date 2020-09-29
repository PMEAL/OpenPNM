function join_by {
    # It was necessary to define this function because
    # Unix's `join` didn't work.
    local IFS="$1"
    shift
    echo "$*"
}

function is_empty {
    if [ "$1" == "- " ]; then
        echo "true"
    else
        echo "false"
    fi
}

function parse_args {
    local temp=()
    for arg in "$@"
    do
        temp+=("\<${arg}\>")
    done
    echo $(join_by "|" "${temp[@]}")
}

function filter_commits_by_label {
    local temp
    local commits=$1    # fetch the first argument
    shift               # removes first arg from list of input args
    temp=$(echo "$commits" | grep -E $(parse_args "$@"))
    temp=$(echo "${temp}" | sed 's/^[ \t]*//; s/[ \t]*$//')
    temp=$(echo "${temp}" | sed -e 's/^/- /')
    echo "$temp"
}

function filter_commits_by_tag_interval {
    local temp
    temp=$(git log "${1}..${2}" -P --author='^((?!GitHub).*)$' --committer=GitHub)
    echo "$temp"
}

function append_to_entry_with_label {
    if [ "$(is_empty "$1")" == "false" ]; then
        echo "### $3"   >> $2
        echo "${1}"     >> $2
        echo ""         >> $2
    fi
}

function get_nth_recent_tag {
    tags=($(git for-each-ref --sort=-creatordate --format '%(refname:strip=2)' refs/tags --count=$1))
    echo "${tags[$(($1-1))]}"
}

# Fetching merge commit messages since the last tag
tag_old=$(get_nth_recent_tag 2)
tag_new=$(get_nth_recent_tag 1)
tag_date=$(git show "$tag_new" --format="%cs")
merge_commits=$(filter_commits_by_tag_interval $tag_old $tag_new)

# Fetching new features/changed API/bugfixes
features=$(filter_commits_by_label "$merge_commits" "Added" "NEW")
enhancements=$(filter_commits_by_label "$merge_commits" "Enhanced" "Optimized" "ENH")
changes=$(filter_commits_by_label "$merge_commits" "Changed" "Removed" "API")
fixes=$(filter_commits_by_label "$merge_commits" "Bugfix" "Hotfix" "Fixed" "BUG")

# Delete "entry" file if already exists
if test -f entry; then
    rm entry
fi

# Delete "CHANGELOG.md" file if already exists
if test -f CHANGELOG.md; then
    rm CHANGELOG.md
fi

# Compile change log
echo -e "## ${tag_new}\n" >> entry
append_to_entry_with_label "$features" entry ":rocket: New features"
append_to_entry_with_label "$enhancements" entry ":cake: Enhancements"
append_to_entry_with_label "$changes" entry ":warning: API changes"
append_to_entry_with_label "$fixes" entry ":bug: Bugfixes"

echo "$(<entry)"

# Modify CHANGELOG.md to reflect new changes
echo -e "# Change log\n" >> CHANGELOG.md
echo "$(<entry)" >> CHANGELOG.md
rm entry

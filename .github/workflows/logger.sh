function join_by {
    # It was necessary to define this function because
    # Unix's `join` didn't work.
    local d=$1
    local f=$2
    shift
    printf %s "$f" "${@/#/$d}"
}

function is_empty {
    if [ "$1" == "- " ]; then
        echo "true"
    else
        echo "false"
    fi
}

function filter_commits_by_label {
    local temp
    local commits=$1            # fetch the first argument
    local exclude_labels=$2     # fetch the second argument
    temp=$(echo "${commits}" | grep -E --ignore-case "$exclude_labels")
    # Strip empty lines (that might include tabs, spaces, etc.)
    temp=$(echo "${temp}" | sed -r '/^\s*$/d')
    # Make each line a bullet point by appending "- " to lines
    temp=$(echo "${temp}" | sed -e 's/^/- /')
    echo "$temp"
}

function filter_commits_exclude_label {
    local temp
    local commits=$1            # fetch the first argument
    local exclude_labels=$2     # fetch the second argument
    # Reverse filter commits by the given labels (i.e. exclude labels)
    temp=$(echo "$commits" | grep -v -E --ignore-case "$exclude_labels")
    # Strip empty lines (that might include tabs, spaces, etc.)
    temp=$(echo "${temp}" | sed -r '/^\s*$/d')
    # Make each line a bullet point by appending "- " to lines
    temp=$(echo "${temp}" | sed -e 's/^/- /')
    echo "$temp"
}

function filter_commits_by_tag_interval {
    local temp
    # --format=%B only outputs commit messages (excluding committer, date, etc.)
    # --first-parent excludes merge commits into the topic branch, ex. dev -> feature
    # temp=$(git log --merges "${1}..${2}" --format=%B --first-parent dev)
    temp=$(git log --merges "${1}..${2}" --format=%B)
    # Remove those merge commits for updating feature branches
    temp=$(echo "${temp}" | grep -v -E "Merge branch")
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

# Fetching new features/enhancements/maintenance/api change/bug fixes/documentation
features=$(filter_commits_by_label "$merge_commits" "#new")
enhancements=$(filter_commits_by_label "$merge_commits" "#enh")
maintenance=$(filter_commits_by_label "$merge_commits" "#maint")
changes=$(filter_commits_by_label "$merge_commits" "#api")
fixes=$(filter_commits_by_label "$merge_commits" "#bug")
documentation=$(filter_commits_by_label "$merge_commits" "#doc")

# Fetching uncategorized merge commits (those w/o keywords)
all_keywords=$(join_by "|" \#new \#enh \#maint \#api \#bug \#doc \# patch \#minor \#major)
uncategorized=$(filter_commits_exclude_label "$merge_commits" "$all_keywords")

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
append_to_entry_with_label "$maintenance" entry ":wrench: Maintenace"
append_to_entry_with_label "$changes" entry ":warning: API changes"
append_to_entry_with_label "$fixes" entry ":bug: Bugfixes"
append_to_entry_with_label "$documentation" entry ":green_book: Documentation"
append_to_entry_with_label "$uncategorized" entry ":question: Uncategorized"

echo "$(<entry)"

# Modify CHANGELOG.md to reflect new changes
echo -e "# Change log\n" >> CHANGELOG.md
echo "$(<entry)" >> CHANGELOG.md
rm entry

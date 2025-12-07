#!/usr/bin/env bash
# shellcheck source=utils.sh
# Generates a changelog from merge commits between the two most recent tags.

set -euo pipefail

log() {
    echo "[INFO] $*" >&2
}

log_section() {
    echo "" >&2
    echo "====== $* ======" >&2
}

join_by() {
    local d="$1"; shift
    local first=1
    for f in "$@"; do
        if (( first )); then
            printf "%s" "$f"
            first=0
        else
            printf "%s%s" "$d" "$f"
        fi
    done
}

is_empty() {
    [[ "${1:-}" == "- " ]] && echo "true" || echo "false"
}

filter_commits_by_label() {
    local commits="$1"
    local label="$2"
    local result
    result=$(echo "$commits" | grep -Ei -- "$label" | \
        sed -E '/^\s*$/d' | \
        sed -E 's/^[[:space:]]*[-]?[[:space:]]*/- /' | \
        sed -E "s/[[:space:]]*$label\b//I" || true)
    local count
    count=$(echo "$result" | grep -c '^- ' || echo "0")
    log "  Found $count commits with label '$label'"
    echo "$result"
}

filter_commits_exclude_label() {
    local commits="$1"
    local exclude_labels="$2"
    local result
    # Filter out labeled commits and common noise patterns from subject lines
    result=$(echo "$commits" | \
        grep -Eiv -- "$exclude_labels" | \
        grep -Eiv -- "^Merge pull request" | \
        sed -E '/^[[:space:]]*$/d' | \
        sed -E 's/^[[:space:]]*[-]?[[:space:]]*/- /' || true)
    local count
    count=$(echo "$result" | grep -c '^- ' || echo "0")
    log "  Found $count uncategorized commits"
    echo "$result"
}

filter_commits_by_tag_interval() {
    local tag_old="$1"
    local tag_new="$2"
    log "Fetching PR merge commits between '$tag_old' and '$tag_new'"
    local result
    # Get full body of merge commits, use separator to split commits
    # Then filter to only PR merges and extract body content (skip merge line)
    result=$(git log --merges "${tag_old}..${tag_new}" --format="%B<<<COMMIT_END>>>" 2>/dev/null | \
        awk '
        BEGIN { RS="<<<COMMIT_END>>>\n?"; ORS="\n" }
        /^Merge pull request/ {
            # Skip the first line (Merge pull request...) and print the rest
            n = split($0, lines, "\n")
            for (i = 2; i <= n; i++) {
                line = lines[i]
                # Skip empty lines, conflict markers, and chore commits
                if (line ~ /^[[:space:]]*$/) continue
                if (line ~ /^#/) continue
                if (line ~ /Conflicts:/) continue
                if (line ~ /^chore\(/) continue
                print line
            }
        }
        ' || true)
    local count
    count=$(echo "$result" | grep -c '.' || echo "0")
    log "Found $count lines from PR descriptions"
    echo "$result"
}

append_to_entry_with_label() {
    local content="$1"
    local file="$2"
    local label="$3"
    if [ "$(is_empty "$content")" = "false" ] && [ -n "$content" ]; then
        log "  Adding section '$label' to $file"
        printf "### %s\n\n%s\n\n" "$label" "$content" >> "$file"
    else
        log "  Skipping empty section '$label'"
    fi
}

# --- Main ---

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/utils.sh"

log_section "Starting changelog generation"

log "Retrieving tags..."
tag_old=$(get_nth_recent_tag 2)
tag_new=$(get_nth_recent_tag 1)
log "Previous tag: $tag_old"
log "Current tag:  $tag_new"

merge_commits=$(filter_commits_by_tag_interval "$tag_old" "$tag_new")

log_section "Categorizing commits by label"
features=$(filter_commits_by_label "$merge_commits" "#new")
enhancements=$(filter_commits_by_label "$merge_commits" "#enh")
maintenance=$(filter_commits_by_label "$merge_commits" "#maint")
changes=$(filter_commits_by_label "$merge_commits" "#api")
fixes=$(filter_commits_by_label "$merge_commits" "#bug")
documentation=$(filter_commits_by_label "$merge_commits" "#docs?")

all_keywords=$(join_by "|" "#new" "#enh" "#maint" "#api" "#bug" "#docs?" "#patch" "#minor" "#major")
uncategorized=$(filter_commits_exclude_label "$merge_commits" "$all_keywords")

log_section "Building changelog entry"

[ -f entry ] && rm -f entry
[ -f CHANGELOG.md ] && rm -f CHANGELOG.md

printf "## %s\n\n" "$tag_new" >> entry
append_to_entry_with_label "$features" entry ":rocket: New features"
append_to_entry_with_label "$enhancements" entry ":cake: Enhancements"
append_to_entry_with_label "$maintenance" entry ":wrench: Maintenance"
append_to_entry_with_label "$changes" entry ":warning: API changes"
append_to_entry_with_label "$fixes" entry ":bug: Bugfixes"
append_to_entry_with_label "$documentation" entry ":green_book: Documentation"
append_to_entry_with_label "$uncategorized" entry ":question: Uncategorized"

log_section "Generated changelog entry"
cat entry

log_section "Writing CHANGELOG.md"
printf "# Change log\n\n" > CHANGELOG.md
cat entry >> CHANGELOG.md
rm -f entry

log "Changelog generation complete"

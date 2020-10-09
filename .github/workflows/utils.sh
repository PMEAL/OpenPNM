function get_nth_recent_tag {
    tags=($(git for-each-ref --sort=-creatordate --format '%(refname:strip=2)' refs/tags --count=$1))
    echo "${tags[$(($1-1))]}"
}


function get_most_recent_tag {
    echo $(git tag | sort -V | tail -1)
}


function get_version {
    version_loc=$1

PYTHON_CODE=$(cat <<EOF
import sys
sys.path.append("$version_loc")
from __version__ import __version__
print(__version__)
EOF
)

    python -c "$PYTHON_CODE"
}

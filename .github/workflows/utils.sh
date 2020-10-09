function get_most_recent_tag {
    echo "$(get_nth_recent_tag 1)"
}


function get_version {
    version_loc=$1
#-- Python code starts here --#
PYTHON_CODE=$(cat <<EOF
import sys
sys.path.append("$version_loc")
from __version__ import __version__
print(__version__)
EOF
)
#--- Python code ends here ---#
    python -c "$PYTHON_CODE"
}


function test_func {
    echo 0123
}


function get_nth_recent_tag {
    git fetch --all --tags
    tags=($(git for-each-ref --sort=-creatordate --format '%(refname:strip=2)' refs/tags --count=$1))
    echo "${tags[$(($1-1))]}"
}

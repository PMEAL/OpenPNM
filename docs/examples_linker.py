import os
from pathlib import Path


def write_nblink(fpath, route):
    rpath = Path(route, fpath)
    new_name = fpath.name.replace(".ipynb", ".nblink")
    fpath = fpath.with_name(new_name)
    with open(fpath, "w") as f:
        f.write("{    \n")
        f.write(f'    "path": "{rpath.as_posix()}"\n')
        f.write("}\n")


dir_docs = Path("./")
dir_examples_root = Path("../")
dir_examples = Path(dir_examples_root, "examples")

for root, subdirectories, files in os.walk(dir_examples):
    for file in files:
        fpath = Path(root, file)
        if ".ipynb_checkpoints" not in str(fpath) and fpath.suffix == ".ipynb":
            # Find path relative to examples root, e.g.: examples/filters/x.ipynb
            rpath = fpath.relative_to(dir_examples_root)
            # Create a similar directory structure under docs/
            rpath.parent.mkdir(parents=True, exist_ok=True)
            # Generate and write .nblink file in the created directory structure
            num_dirs_deep = len(rpath.parent.parts)
            route = Path("../" * num_dirs_deep, dir_examples_root)
            write_nblink(rpath, route)

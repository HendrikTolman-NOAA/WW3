import os
import glob
import shutil

def modify_compilers():
    # Find full path to icpx if available, else fallback to 'icpx'
    icpx_path = shutil.which('icpx') or 'icpx'

    # Find all compilers.yaml files and the active environment's spack.yaml file
    search_paths = [
        os.path.expanduser('~/.spack/linux/compilers.yaml'),
        os.path.expanduser('~/.spack/*/compilers.yaml'),
        'spack/var/spack/environments/ww3-intel/spack.yaml',
        '../spack/var/spack/environments/ww3-intel/spack.yaml',
        '/home/runner/work/WW3/WW3/spack/var/spack/environments/ww3-intel/spack.yaml'
    ]

    paths = []
    for pattern in search_paths:
        paths.extend(glob.glob(pattern))

    for p in set(paths):
        if not os.path.exists(p):
            continue
        print(f"Checking compiler configuration file: {p}")
        with open(p, 'r') as f:
            lines = f.readlines()

        new_lines = []
        in_intel = False
        modified = False

        for line in lines:
            # Check if we are entering an intel or oneapi compiler block
            if 'spec: intel' in line or 'spec: oneapi' in line:
                in_intel = True
            elif 'spec:' in line or '- compiler:' in line:
                # If we see another compiler section, reset the flag unless it's intel-specific
                if 'intel' not in line and 'oneapi' not in line:
                    in_intel = False

            # If we are in an intel/oneapi compiler block, modify cxx path to icpx
            if in_intel and 'cxx:' in line:
                indent = len(line) - len(line.lstrip())
                new_line = ' ' * indent + f'cxx: {icpx_path}\n'
                if line != new_line:
                    print(f"Modifying cxx path in {p}:")
                    print(f"  Old: {line.strip()}")
                    print(f"  New: {new_line.strip()}")
                    line = new_line
                    modified = True
            new_lines.append(line)

        if modified:
            with open(p, 'w') as f:
                f.writelines(new_lines)
            print(f"Successfully updated: {p}")

if __name__ == '__main__':
    modify_compilers()

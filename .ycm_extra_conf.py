import os
import re
import ast

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))

flags = []
pattern = re.compile('^([A-Z0-9]+)\s*=\s*(\[.*\])')
with open(os.path.join(SCRIPT_PATH, "config.py")) as f:
    for line in f:
        result = pattern.search(line)
        if result is None:
            continue
        key = result.group(1)
        if 'LINK' in flags or (not 'FLAGS' in key and not 'INCLUDE' in key):
            continue
        try:
            result = ast.literal_eval(result.group(2))
            if 'INCLUDE' in key:
                flags += [ '-I' + path.replace('#', '') for path in result ]
            else:
                flags += result
        except:
            pass
        
def MakeRelativePathsInFlagsAbsolute(flags, working_directory):
    if not working_directory:
        return list(flags)
    new_flags = []
    make_next_absolute = False
    path_flags = ['-isystem', '-I', '-iquote', '--sysroot=']
    for flag in flags:
        new_flag = flag

        if make_next_absolute:
            make_next_absolute = False
            if not flag.startswith('/'):
                new_flag = os.path.join(working_directory, flag)

        for path_flag in path_flags:
            if flag == path_flag:
                make_next_absolute = True
                break

            if flag.startswith(path_flag):
                path = flag[len(path_flag):]
                new_flag = path_flag + os.path.join(working_directory, path)
                break

        if new_flag:
            new_flags.append(new_flag)
    return new_flags

def FlagsForFile(filename, **kwargs):
    final_flags = MakeRelativePathsInFlagsAbsolute(flags, SCRIPT_PATH) + ['-x', 'c++']
    return {
        'flags': final_flags,
        'do_cache': True
    }

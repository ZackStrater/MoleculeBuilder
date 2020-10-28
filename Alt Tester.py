from termcolor import cprint
import itertools

a = {"a": 1, "b": 2, "c": 3}
b = {"d": 4, "e": 5, "f": 6}
c = {"g": 7, "h": 8, "i": 9}

def hierarchy_check(*libraries):
    def check_func(index, library):
        truncated_lib = itertools.islice(library.items(), index)
        for k, v in truncated_lib:
            print(k, v)
    total_library = {}
    for lib in libraries:
        total_library.update(lib)
    for i, (key, value) in enumerate(total_library.items()):
        cprint("new entry:", "blue")
        print(i, key, value)
        cprint("checked entries:", "blue")
        check_func(i, total_library)
        print("\n")

hierarchy_check(a, b, c)


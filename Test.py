
test = "C12=CC=CN=C1{cccc}"

import re
def replace(matchobj):
    first_edit = re.sub(r"([a-z])", r"\1*", matchobj.group())
    return re.sub(r"[{}]", r"", first_edit)

print(re.sub(r"({)([a-z]+)(})", replace, test))


test2 = "C12=CC=CN=C1C=CC=C2}"
if re.search(r"[{}]", test2):
    print("yes")

from collections import Counter

x = {1, 2, 3}

y = {1: 2, 2: 2, 3: 5, 4: 8, 5: 7}

for num in x:
    y.pop(num, None)

print(y)
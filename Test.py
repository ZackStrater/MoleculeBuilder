import random

def recursive_func(list):
    list.append(random.randint(3,20))
    print(f"this is the list: {list}")
    if len(list) > 20:
        return
    recursive_func(list)


recursive_func([3])
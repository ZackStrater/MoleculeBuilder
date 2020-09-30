
x = [[1, 1, 11], [1, 1, 1], [0, 2, 2]]


def add_list(list):
    def double_odd_list(list):
        if sum(list) % 2 != 0:
            return double_odd_list([x*2 for x in list])
        else:
            return sum(list)
    return double_odd_list(list)


x.sort(key=add_list, reverse=True)

print(x)
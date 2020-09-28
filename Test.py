
x = [1, 2, 3, 4, 5, 6, 7, 8]


for num in x:
    if num%2 == 0:
        x.remove(num)
    elif num%3 == 0:
        x.remove(num)
print(x)
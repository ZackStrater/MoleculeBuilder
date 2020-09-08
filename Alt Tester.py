from collections import Counter

listA = [1, 2, 2, 2, 3, 7]

listB = [1, 2, 6, 3, 4, 2, 2]


def listA_in_listB(listA, listB):

    cA = Counter(listA)
    cB = Counter(listB)

    for key in cA:
        if key not in cB or cA[key] > cB[key]:
            print("nope")
            return

    print("ListA in listB")

listA_in_listB(listA, listB)
a = [1,2,3,4,5]
b = [2,3,6,9]
def returnNotMatches(a, b):
    return [[x for x in a if x not in b], [x for x in b if x not in a]]


print(returnNotMatches(a,b))

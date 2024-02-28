def ordset(n,x=0):
    for _ in range(n+1):
        yield x
        x = (x or[]) + [x]

number = 3
a = ordset(number, x=[])

i=0
print('[', end='')
while i < number:
    print(next(a), end = '')
    if not i == number-1:
        print(',', end = '')
    i+=1
print(']')
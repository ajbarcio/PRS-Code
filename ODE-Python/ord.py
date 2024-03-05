import time
def ordset(n,x=0):
    for _ in range(n+1):
        yield x
        x = (x or[]) + [x]
        # time.sleep(5.5)

number = 7

a = ordset(number, x=[])

i=0
print('[', end='')
while i < number:
    print(next(a), end = '')
    if not i == number-1:
        print(',', end = '')
    i+=1
print(']')

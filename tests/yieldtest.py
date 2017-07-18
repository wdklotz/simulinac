def genfunc():
    print('in body of genfunc')
    for k in range(10):
        print('loop in genfunc\t k=',k)
        yield k

for i in genfunc():
    print('loop in main \t i=',i)

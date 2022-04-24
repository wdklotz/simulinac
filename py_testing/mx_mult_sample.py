import numpy as np
v=np.array([1,1])
print ('v >>\n',v)

a=np.array([[1,-1],[1,1]])
print('a >>\n',a)
av= a.dot(v)
print('a*v=a.dot(v) >>\n',av)

b=np.array([[1,2],[3,4]])
print('b >>\n',b)
bv= b.dot(v)
print('b*v=b.dot(v) >>\n',bv)

print ('a*b=a.dot(b) >>\n',a.dot(b))
p=np.einsum('ij,jk',a,b)
print ('np.einsum(\'ij,jk\',a,b) >>\n',p)
abv= p.dot(v)
print ('(a*b)*v >>\n',abv)

print('b*a=b.dot(a) >>\n',b.dot(a))
p=np.einsum('ij,jk',b,a)
print('np.einsum(\'ij,jk\',b,a) >>\n',p)
bav= p.dot(v)
print('(b*a*)v >>\n',bav)



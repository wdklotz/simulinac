import h5py
import numpy as np

arr=np.zeros((4,6))
for i in range(4):
    for j in range(6):
        arr[i][j]=i*6+j+1

print('open "file.h5"')
file = h5py.File('file.h5','w')
print('create dataset "/dset"')
dataset = file.create_dataset('dset',(4,6),data=arr)
data_read=dataset[...]
file.close()
print('close "file.h5"')
print(data_read,'\n')

print('open "file.h5"')
file = h5py.File('file.h5','r+')
print('create dataset "/Mygroupt/dest"')
group = file.create_group('MyGroup')
dataset1 = file.create_dataset('/MyGroup/dset',(4,6),data=2*arr)
data_read1=dataset1[...]
file.close()
print('close "file.h5"')
print(data_read1,'\n')

print('open "file.h5"')
file = h5py.File('file.h5','r')
print('read dataset "/dest"')
dataset=file['/dset']
data_read=dataset[...]
print(data_read,'\n')

print('read dataset "/MyGroup/dest"')
dataset1=file['/MyGroup/dset']
data_read1=dataset1[...]
print(data_read1)
print('close "file.h5"')
file.close()

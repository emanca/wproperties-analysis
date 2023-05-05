import os

for i in range(100):
    print('running iteration', i)
    os.system('python boost2pandas.py -iteration {}'.format(i)) 

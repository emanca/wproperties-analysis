import os

for i in range(100):
    print('\n\n\n\n\n\n\n\\n\n\n')
    print('running iteration', i)
    os.system('python boost2pandas_twocharges.py -iteration {}'.format(i))
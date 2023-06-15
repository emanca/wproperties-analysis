import os

for i in range(100):
    print('\n\n\n\n\n\n\n\n\n\n\n\n\n')
    print('running fits, iteration',i)
    os.system('combinetf.py --fitverbose 9 -t0 --seed 260292 --yieldProtectionCutoff 100. --allowNegativePOI --doh5Output ../templateMaker/booststrap_tests_redstat/Wlike_iteration_{}_twocharges_redstat.hdf5 -o FitRes/bootstrap_tests_redstat/fit_Wlike_iteration_{}_twocharges_redstat.root --preconditioning'.format(i,i))
for j in range(100):
    print('\n\n\n\n\n\n\n\n\n\n\n\n\n')
    print('running fits, iteration',j)
    os.system('combinetf.py --fitverbose 9 -t0 --seed 260292 --yieldProtectionCutoff 100. --allowNegativePOI --doh5Output ../templateMaker/bootstrap_tests_redstat/Wlike_iteration_{}_twocharges_redstat.hdf5 -o FitRes/bootstrap_tests_redstat/fit_Wlike_iteration_{}_twocharges_redstatBBB.root --preconditioning -- binByBinStat'.format(j,j))

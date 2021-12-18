import ROOT
import math
import os
import sys
import copy

#  /scratch/sroychow/mcsamples_2016-V1MCFinal.txt
# /QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 3616000.0
# /QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 3160000.0
# /QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 1655000
# /QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 448700
# /QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 105200
# /QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 105200
# /QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 25510
# /QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 25510
# /QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 8624
# /QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 8624
# /QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 8624
# /QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 792.2
# /QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 792.2
# /QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM, 792.2
# /QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 79.18
# /QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 79.18
# /QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM, 79.18
# /QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 25.22
# /QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 25.22
# /QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 25.22
# /QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 4.717
# /QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 4.717
# /QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM, 4.717
# /QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM, 1.615
# /QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM, 1.615



# /QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 3616000.0
# /QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 3160000.0
# /QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 1655000
# /QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 448700
# /QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 105200
# /QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 105200
# /QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 25510
# /QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 8624
# /QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 8624
# /QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 792.2
# /QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 792.2
# /QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext2-v1, 792.2
# /QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 79.18
# /QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 79.18
# /QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext2-v1, 79.18
# /QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 25.22
# /QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 25.22
# /QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 4.717
# /QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 4.717
# /QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext2-v1, 4.717
# /QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/-v1, 1.615
# /QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/_ext1-v1, 1.615




sampleList = []
sampleList.append("QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
sampleList.append("QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
sampleList.append("QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
sampleList.append("QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
sampleList.append("QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
sampleList.append("QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM")
sampleList.append("QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
# sampleList.append("QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM")
sampleList.append("QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM") 
# sampleList.append("QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM") 
# sampleList.append("QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM") 
sampleList.append("QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM") 

path = '/gpfs/ddn/srm/cms/store/user/sroychow/NanoAOD2016-V2/'
# /gpfs/ddn/srm/cms/store/user/sroychow/NanoAOD2016-V2/QCD_Pt*/*/*/*/*.root
for s in sampleList :
    sName = s.split('/')[0]
    sVerFull =s.split('/')[1]
    sVer = sVerFull.replace('RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7','')
    sVer = sVer.replace('RunIISummer16NanoAODv6-PUMoriond17_Nano25Oct2019_backup_102X_mcRun2_asymptotic_v7','')
    sVer = sVer.replace('-v1','')
    # print(sName+sVer)
    
    os.system('mkdir -p '+sName+sVer)
    outName = sName+sVer+'/tree.root'
    # print('hadd -f '+ outName +' '+path+sName+'/'+sVerFull+'/*/*/*.root')
    os.system('hadd -f '+ outName +' '+path+sName+'/'+sVerFull+'/*/*/*.root')
    
    # pathNT = '/scratchssd/bertacch/wmass/CMSSW_10_2_9/src/PhysicsTools/NanoAODTools/scripts/'
    # os.system(pathNT+'haddnano.py '+ outName +' '+path+sName+'/'+sVerFull+'/*/*/*.root')

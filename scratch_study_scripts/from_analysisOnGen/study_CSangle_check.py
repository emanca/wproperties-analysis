import ROOT 
import math 
import copy

def azimuth(phi):	
    if phi<0.0:
        phi += 2*math.pi	
    return phi

def evalCSangle_Lorenzo(muon,neutrino) :

    Lp4  = ROOT.TLorentzVector(0.,0.,0.,0.)
    Np4 = ROOT.TLorentzVector(0.,0.,0.,0.)
    # if hasattr(muon, "pt"):
    #     Lp4.SetPtEtaPhiM(muon.pt, muon.eta, muon.phi, muon.mass)
    #     Np4.SetPtEtaPhiM(neutrino.pt, neutrino.eta, neutrino.phi, 0.)
    # else:
    Lp4.SetPtEtaPhiM(muon[0], muon[1], muon[2], muon[3])
    Np4.SetPtEtaPhiM(neutrino[0], neutrino[1], neutrino[2], 0.)
    Wp4 = Lp4 + Np4
    
    Wp4_rot = copy.deepcopy(Wp4)
    Lp4_rot = copy.deepcopy(Lp4)

    # align W/L along x axis
    Wp4_rot.RotateZ( -Wp4.Phi() )
    Lp4_rot.RotateZ( -Wp4.Phi() )

    # first boost
    boostL = Wp4_rot.BoostVector()
    boostL.SetX(0.0)
    boostL.SetY(0.0)
    Lp4_rot.Boost( -boostL )
    Wp4_rot.Boost( -boostL )

    # second boost
    boostT = Wp4_rot.BoostVector()
    Lp4_rot.Boost( -boostT )

    # the CS frame defines the z-axis according to the W pz in the lab 
    flip_z = -1 if Wp4.Rapidity()<0.0 else +1
    # flip_z=1
    
    phi = math.atan2(Lp4_rot.Py()*flip_z, Lp4_rot.Px())
    # if phi<0: phi = phi + 2*math.pi


    # compute PS point
    #print 'flip', flip_z, ': ', azimuth(Lp4_rot.Phi()*flip_z), ' --- ', math.atan2(Lp4_rot.Py()*flip_z,Lp4_rot.Px()*flip_z)
    # return [Lp4_rot.CosTheta()*flip_z, math.atan2(Lp4_rot.Py()*flip_z, Lp4_rot.Px())] #official
    return [Lp4_rot.CosTheta()*flip_z, phi] #official
    
    
    
    # return [Lp4_rot.CosTheta()*flip_z, math.atan2(Lp4_rot.Py()*flip_z, Lp4_rot.Px()*flip_z)]
    # return [Lp4_rot.CosTheta()*flip_z, azimuth(Lp4_rot.Phi()*flip_z)]	
    
def fixPhi(phi,y) :
    if y>=0. : 
        out = phi
    if phi>=0. : 
        out = math.pi-phi
    else : 
        out = -math.pi-phi
    return out 

def evalCSangle(muon, neutrino):
    	
    m = ROOT.TLorentzVector()
    n = ROOT.TLorentzVector()
    w = ROOT.TLorentzVector()
    
    m.SetPtEtaPhiM(muon[0], muon[1], muon[2], 0.105)
    n.SetPtEtaPhiM(neutrino[0],neutrino[1],neutrino[2], 0.)
    
    w = m + n

    sign  = abs(w.Z())/w.Z()
    
    ProtonMass = 0.938
    BeamEnergy = 6500.000
    
    p1 = ROOT.TLorentzVector()
    p2 = ROOT.TLorentzVector()
    
    p1.SetPxPyPzE(0, 0, sign*BeamEnergy, math.sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass)) 
    p2.SetPxPyPzE(0, 0, -1*sign*BeamEnergy, math.sqrt(BeamEnergy*BeamEnergy+ProtonMass*ProtonMass))
    
    print("pre, p1 x,y=", p1.X(), p1.Y())
    wp = ROOT.TLorentzVector()
    wp.SetPxPyPzE(0.5, 0., 0., 10)
    # print("w", w.X(), w.Y(), w.Z(),w.E()) 
    # print(wp.BoostVector().X(),wp.BoostVector().Y(),wp.BoostVector().Z())
    p1.Boost(-wp.BoostVector())
    p2.Boost(-wp.BoostVector())
    # p1.Boost(-w.BoostVector())
    # p2.Boost(-w.BoostVector())
    print("post, p1 x,y=", p1.X(), p1.Y())
    print("post, p2 x,y=", p2.X(), p2.Y())
    
    CSAxis = (p1.Vect().Unit()-p2.Vect().Unit()).Unit() #quantise along axis that bisects the boosted beams
    
    yAxis = (p1.Vect().Unit()).Cross((-p2.Vect().Unit())) #other axes #FIX ADDED, -p2.   
    yAxis = yAxis.Unit()
    xAxis = yAxis.Cross(CSAxis)
    xAxis = xAxis.Unit()
    print("xax", xAxis.X(), xAxis.Y(), xAxis.Z(),)

    m.Boost(-w.BoostVector())

    phi = math.atan2((m.Vect()*yAxis),(m.Vect()*xAxis))
    if phi<0: phi = phi + 2*math.pi
    
    return math.cos(m.Angle(CSAxis)), phi
    


inFile_raw = ROOT.TFile.Open("/scratchnvme/wmass/WJetsNoCUT_v2/tree_0_1.root")
# inFile_raw = ROOT.TFile.Open("/scratchnvme/wmass/NanoAOD2016-V2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext2/tree.root")
tree = inFile_raw.Get("Events")
Nmax=100

for i in range(1,Nmax) :
    mu = []
    nu = []
    tree.GetEntry(i)
    
    muInd = tree.GenPart_preFSRMuonIdx
    nuInd = tree.GenPart_NeutrinoIdx
    
    if muInd==-99 : continue
    
    mu.append(tree.GenPart_pt[muInd])
    mu.append(tree.GenPart_eta[muInd])
    mu.append(tree.GenPart_phi[muInd])
    mu.append(tree.GenPart_mass[muInd])
    
    nu.append(tree.GenPart_pt[nuInd])
    nu.append(tree.GenPart_eta[nuInd])
    nu.append(tree.GenPart_phi[nuInd])
    
    CSCosTheta_eval, CSphi_eval = evalCSangle_Lorenzo(mu,nu)
    CSCosTheta_eval_eli, CSphi_eval_eli = evalCSangle(mu,nu)
    
    CSCosTheta_tree = tree.CStheta_preFSR
    CSphi_tree = tree.CSphi_preFSR-math.pi
    yW = tree.Wrap_preFSR
    
    CSCosTheta_tree = CSCosTheta_eval_eli
    CSphi_tree = CSphi_eval_eli
    
    
    # CSphi_eval = fixPhi(CSphi_eval,yW)
    # CSphi_tree = fixPhi(CSphi_tree,yW)

    # print("eval-tree (costheta, phi)=", CSCosTheta_eval-CSCosTheta_tree, CSphi_eval-CSphi_tree, " ctheta=",CSCosTheta_eval, " phi=",CSphi_eval, " ctheta_t=", CSCosTheta_tree, " cphi_t=",CSphi_tree)
    
    
    
    
    




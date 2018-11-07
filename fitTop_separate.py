#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import *
import tdrstyle
import numpy as np 
tdrstyle.setTDRStyle()

fVar="Puppijet0_msd"
fPt="Puppijet0_pt"
fCut="(Puppijet0_pt>350)"
fWeight="weight"
fLumi=41.1
fBinWidth =5
fDir="/afs/cern.ch/user/c/cmantill/work/public/forJeff/sklimWtag"
fPlotDir="plots/"
fOutput="top.root"

fFiles = {'data': 'Data',
          #'wlnu_mc' : 'WJets',
          #'wlnu_data': 'WJets',
          'st_mc':   'ST',
          'st_data':   'ST',
          'tt_mc':   'TT',
          'tt_data':   'TT',
          'tt_fakeW_mc': 'TT',
          'tt_realW_mc': 'TT',
          'tt_fakeW_data': 'TT',
          'tt_realW_data': 'TT',
          'mc': 'PseudoData'
          }

fPDFs = {
         #'wlnu_mc':       'Exp_mc',
         #'wlnu_data':     'Exp_data',
         #'st_mc':        'GausErfExp_tt',
         'st_mc':         'PolyGaus_st_mc',
         'st_data':       'PolyGaus_st_data',
         'tt':            'GausErfExp_tt',
         'tt_mc':         'GausErfExp_tt_mc',
         'tt_data':       'GausErfExp_tt_data',
         #'tt_data':       'GausGaus_tt_data',
         'tt_fakeW_mc':   'GausErfExp_tt_mc',
         'tt_realW_mc':   'GausErfExp_tt_mc',
         'tt_fakeW_data': 'GausErfExp_tt_data',
         'tt_realW_data': 'GausErfExp_tt_data',
         }

r.gSystem.Load("./PDFs/HWWLVJRooPdfs_cxx.so")
r.gSystem.Load("./PDFs/PdfDiagonalizer_cc.so")
r.gSystem.Load("./PlotStyle/Util_cxx.so")
r.gSystem.Load("./PlotStyle/PlotUtils_cxx.so")

def parser():
    parser = OptionParser()
    parser.add_option('--tag',     action='store', type='string', dest='tag',   default='CA15TopSelTau32', help='samples tag')
    parser.add_option('--xMin',    action='store', type='float',  dest='xmin',  default=40,                help='x-min')
    parser.add_option('--xMax',    action='store', type='float',  dest='xmax',  default=250,               help='x-max')
    parser.add_option('-b',        action='store_true',           dest='noX',   default=True,              help='no X11 windows')
    parser.add_option('--match',   action='store_true',           dest='match', default=False,             help='divide tt sample in matched and unmatched W')
    (options,args) = parser.parse_args()
    return options

# can do a couple of things here: colors, errors, add legends!
def drawFrame(iFrame,iData,iFuncs,isPoisson=False):
    if isPoisson:
        iData.plotOn(iFrame,r.RooFit.DataError(r.RooAbsData.Poisson),r.RooFit.XErrorSize(0))
    else:
        iData.plotOn(iFrame)
    iColor=50
    for pFunc in iFuncs:
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1))

def draw(iVar,iData,iFuncs,iRooFitResult,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    drawFrame(lFrame,iData,iFuncs)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %f GeV"%fBinWidth);
    lchisqndof = getChi2NDOF(iData,iFuncs,iRooFitResult,iVar)
    addInfo = getPavetext()
    addInfo.AddText("#chi^{2}/nDOF = %.3f"%lchisqndof)
    print 'Chi^2/NDOF %.3f'%lchisqndof
    lFrame.Draw()
    addInfo.Draw()
    lCan.Modified()
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")

def getPavetext():
    addInfo = TPaveText(0.2,0.75,0.5,0.9,"NDC")
    addInfo.SetFillColor(1)
    addInfo.SetLineColor(1)
    addInfo.SetFillStyle(1)
    addInfo.SetBorderSize(0)
    addInfo.SetTextFont(42)
    addInfo.SetTextSize(0.040)
    addInfo.SetTextAlign(12)
    return addInfo

def drawDataMc(iVar,iData,iFuncsData,iMc,iFuncsMc,iRooFitResult_mc,iRooFitResult_data,params,iLabel='A'):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    drawFrame(lFrame,iData,iFuncsData,True)
    #iData.plotOn(lFrame)
    #iMc.plotOn(lFrame)
    #iFuncsData[0].plotOn(lFrame,r.RooFit.LineColor(1))
    #iFuncsMc[0].plotOn(lFrame,r.RooFit.LineColor(2))
    data_chisqndof = getChi2NDOF(iData,iFuncsData,iRooFitResult_mc,iVar)
    mc_chisqndof   = getChi2NDOF(iMc,iFuncsMc,iRooFitResult_data,iVar)

    print 'scale factor + resolution for t sample: '
    print 'scale =', params["data_mean_val"]/params["mc_mean_val"], "+/-", params["data_mean_val"]/params["mc_mean_val"]*np.sqrt((params["mc_mean_err"]/params["data_mean_val"])**2 + (params["mc_mean_err"]/params["mc_mean_val"])**2)
    print 'resolution =', params["data_sigma_val"]/params["mc_sigma_val"], "+/-", params["data_sigma_val"]/params["mc_sigma_val"]*np.sqrt((params["data_sigma_err"]/params["data_sigma_val"])**2 + (params["mc_sigma_err"]/params["mc_sigma_val"])**2)
    print  
    print 'chisq/ndof for data =', data_chisqndof
    print 'chisq/ndof for mc =', mc_chisqndof
    drawFrame(lFrame,iMc,iFuncsMc)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %f GeV"%fBinWidth);

    addInfo = getPavetext()
    addInfo.AddText("Data #chi^{2}/nDOF = %.3f"%data_chisqndof) 
    addInfo.AddText("MC #chi^{2}/nDOF = %.3f"%mc_chisqndof) 
    lFrame.Draw()
    addInfo.Draw() 
    lCan.Modified() 
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")

def getChi2NDOF(iData,iModel,iRooFitResult,iVar):
    lDataHist = iData.binnedClone(iData.GetName()+"_binnedClone",iData.GetName()+"_binnedClone");
    pChi2     = iModel[0].createChi2(lDataHist,r.RooFit.Range(100.,250.),r.RooFit.Extended(r.kTRUE),r.RooFit.DataError(r.RooAbsData.Poisson))
    #print '#chi sq',pChi2.getVal()
    #print 'number of bins', iVar.getBins()
    #print 'number of free parameters', iRooFitResult.floatParsFinal().getSize()
    return pChi2.getVal()/(int(iVar.getBins() - iRooFitResult.floatParsFinal().getSize())); 

def getPull(iVar,iPlot):
    lHPull = iPlot.pullHist();
    pPull  = iVar.frame(r.RooFit.Title("Pull"), r.RooFit.Bins(int(iVar.getBins())));
    lLine  = r.TLine(iVar.getMin(),0,iVar.getMax(),0);
    lLine.SetLineWidth(2); lLine.SetLineColor(r.kRed);
    pPull.addObject(lLine);
    pPull.addPlotable(lHPull,"P");
    pPull.SetTitle("");
    pPull.GetXaxis().SetTitle("");
    pPull.GetYaxis().SetTitle("pull");
    pPull.GetYaxis().SetRangeUser(-4,4);
    pPull.GetXaxis().SetTitleSize(0.10);
    pPull.GetXaxis().SetLabelSize(0.10);
    pPull.GetYaxis().SetTitleSize(0.10);
    pPull.GetYaxis().SetLabelSize(0.10);
    return pPull

# add gaussian constraint                                                                                                   
def addConstraint(iWorkspace,iVar,iMean,iSigma,iList):
    print '---- Adding gaussian constraint to %s'%iVar.GetName()
    lMean = r.RooRealVar("%s_mean"%iVar.GetName(),"%s_mean"%iVar.GetName(),iMean);
    lSigma = r.RooRealVar("%s_sigma"%iVar.GetName(),"%s_sigma"%iVar.GetName(),iSigma);
    lConstraint = r.RooGaussian("constraint_pdf_%s"%iVar.GetName(),"constraint_pdf_%s"%iVar.GetName(),iVar,lMean,lSigma)
    lConstraint.Print()
    iList.append(lConstraint.GetName())
    getattr(iWorkspace,"import")(lConstraint,r.RooFit.RecycleConflictNodes())


# make Pdf from model, probably should put parameters in dictionary
def makePdf(iWorkspace,iLabel,iModel,iMc=False):
    print '---- Making pdf' 
    lVar = iWorkspace.var(fVar);
    lModelPdf = None
    lTag = "%s_%s"%(iModel,iLabel)
    pVarHigh       = r.RooRealVar("%sHi_%s"%(lVar.GetName(),lTag),"%sHi_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    pVarHigh1      = r.RooRealVar("%sHi1_%s"%(lVar.GetName(),lTag),"%sHi1_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    lConstraints = []


    if iModel == "Exp_mc":
        lC_Exp_mc        = r.RooRealVar("c_mc"+lTag,"c_mc"+lTag,-0.03,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_mc);
        
    if iModel == "Exp_data":
        lC_Exp_data      = r.RooRealVar("c_data"+lTag,"c_data"+lTag,-0.03,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_data);
        
    if "PolyGaus_st" in iModel:
        lMean_Gaus       = r.RooRealVar("mean"+lTag,"mean"+lTag,175.,120.,190.)
        lSigma_Gaus      = r.RooRealVar("sigma"+lTag,"sigma"+lTag,15.,0.,100.)
        pGaus            = r.RooGaussian("gaus"+lTag,"gaus"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        lMean2_Gaus      = r.RooRealVar("mean2"+lTag,"mean2"+lTag,150.,140.,175.)
        lSigma2_Gaus     = r.RooRealVar("sigma2"+lTag,"sigma2"+lTag,80.,40.,100.)
        pGaus2           = r.RooGaussian("gaus2"+lTag,"gaus2"+lTag,lVar,lMean2_Gaus,lSigma2_Gaus)
        la1              = r.RooRealVar("a1"+lTag,"a1"+lTag,-0.003,-0.005,0.)
        #la1             = r.RooRealVar("a1"+lTag,"a1"+lTag,-0.04,-.1,0.)
        #la2             = r.RooRealVar("a2"+lTag,"a2"+lTag,0.1,0.,1.)
        #pPoly           = r.RooExponential("poly"+lTag,"poly"+lTag,lVar,la1)
        pPoly            = r.RooPolynomial("poly"+lTag,"poly"+lTag,lVar,r.RooArgList(la1),1)
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pPoly),r.RooArgList(pVarHigh))
        if 'data' in iModel: # or 'mc' in iModel:
            addConstraint(iWorkspace,lMean_Gaus,177,5,lConstraints)
            addConstraint(iWorkspace,lSigma_Gaus,17,5,lConstraints)
            addConstraint(iWorkspace,la1,-0.002,0.0008,lConstraints)

    if "ErfExpGaus_st" in iModel:
        lC_ErfExp        = r.RooRealVar("c"+lTag,"c"+lTag,-0.05,-0.2,0.0)
        lOffSet_ErfExp   = r.RooRealVar("offset"+lTag,"offset"+lTag,175.,140.,300.)
        lWidth_ErfExp    = r.RooRealVar("width"+lTag,"width"+lTag,30.,10.,300.)
        lMean_Gaus       = r.RooRealVar("mean"+lTag,"mean"+lTag,175.,140.,200.)
        lSigma_Gaus      = r.RooRealVar("sigma"+lTag,"sigma"+lTag,7.,0.,40.)
        pGaus            = r.RooGaussian("gaus"+lTag,"gaus"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        #add W 
        #lMean_Gaus2     = r.RooRealVar("mean2"+lTag,"mean2"+lTag,85.,75.,95.)
        #lSigma_Gaus2    = r.RooRealVar("sigma2"+lTag,"sigma2"+lTag,10.,2.,20.)
        #pGaus2          = r.RooGaussian("gaus2"+lTag,"gaus2"+lTag,lVar,lMean_Gaus2,lSigma_Gaus2)
        
        pErfExp          = r.RooErfExpPdf("erfExp"+lTag,"erfExp"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        #lModelPdf       = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus,pGaus2),r.RooArgList(pVarHigh,pVarHigh1));
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus),r.RooArgList(pVarHigh));

    if iModel == "ErfExp_st_mc":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if iModel == "ErfExp_st_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if  "GausErfExp_tt" in iModel:
        print iLabel
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,0.,200.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,30.,0.,200.)
        lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,130.,200.)
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        if iLabel == "tt_mc" or iLabel == "tt_data":
            lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,20.,0.,200.)
        if iLabel == "tt_mc":
            lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,10.,5.,15.)
        # if "realW" in iLabel:
        #     lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,140.,200.)
        #     lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,15.,0.,30.)
        #     lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,250.,0.,1000.)
        #     lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,60.,0.,200.)
        if "fakeW" in iLabel:  
            lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,100.,-1000.,1000.)
            lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,30.,-1000.,1000.)
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,125.,100.,200.)
            lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        pVarFrac       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        #lC_ErfExp.setConstant(r.kTRUE);
        #lOffSet_ErfExp.setConstant(r.kTRUE);
        #lWidth_ErfExp.setConstant(r.kTRUE);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp),r.RooArgList(pVarHigh),1);
        if 'mc' in iLabel and not 'fakeW' in iLabel and iLabel != 'tt_mc':
            addConstraint(iWorkspace,lC_ErfExp,177,5,lConstraints)
            addConstraint(iWorkspace,lMean_Gaus,177,10,lConstraints)
            addConstraint(iWorkspace,lSigma_Gaus,15,10,lConstraints)
            #addConstraint(iWorkspace,lWidth_ErfExp,60,20,lConstraints)
            addConstraint(iWorkspace,lOffSet_ErfExp,200,10,lConstraints)
            #Puppijet0_msdHi_GausErfExp_tt_data_tt_realW_data    5.8297e-01 +/-  1.75e-01
            #c_GausErfExp_tt_data_tt_realW_data   -2.4649e-02 +/-  3.15e-02
            #mean_GausErfExp_tt_data_tt_realW_data    1.7787e+02 +/-  3.53e-01
            #number_tt_realW_data    8.2941e+03 +/-  6.27e+01
            #offset_GausErfExp_tt_data_tt_realW_data    2.0000e+02 +/-  1.71e-01
            #sigma_GausErfExp_tt_data_tt_realW_data    1.5277e+01 +/-  2.61e+00
            #width_GausErfExp_tt_data_tt_realW_data    6.1304e+01 +/-  2.23e+01
            addConstraint(iWorkspace,lOffSet_ErfExp,-4.2016e-02,0.01,lConstraints)
        if 'data' in iLabel and not 'fakeW' in iLabel: # and not 'realW' in iLabel:
            #addConstraint(iWorkspace,lMean_Gaus,175,10,lConstraints)
            #addConstraint(iWorkspace,lSigma_Gaus,15,10,lConstraints)
            #addConstraint(iWorkspace,lWidth_ErfExp,10,5,lConstraints)
            addConstraint(iWorkspace,lOffSet_ErfExp,200,10,lConstraints)
            addConstraint(iWorkspace,lOffSet_ErfExp,-4.2016e-02,0.01,lConstraints)

        #if iLabel == 'tt_mc':
        #    addConstraint(iWorkspace,lMean_Gaus,173,10,lConstraints)
        #    addConstraint(iWorkspace,lSigma_Gaus,11,6,lConstraints)
        #    addConstraint(iWorkspace,lOffSet_ErfExp,200,10,lConstraints)
        if iLabel == 'tt_data':
            addConstraint(iWorkspace,lMean_Gaus,173,10,lConstraints)
            addConstraint(iWorkspace,lSigma_Gaus,11,6,lConstraints)
            addConstraint(iWorkspace,lOffSet_ErfExp,200,10,lConstraints)

  # Puppijet0_msdHi_GausErfExp_tt_mc_tt_mc    2.5475e-01 +/-  1.74e-01
  # Puppijet0_msdHi_PolyGaus_st_mc_st_mc    8.2245e-01 +/-  5.56e-01
  # a1PolyGaus_st_mc_st_mc   -1.0582e-06 +/-  3.98e-02
  # c_GausErfExp_tt_mc_tt_mc   -2.5316e-02 +/-  8.73e-03
  # meanPolyGaus_st_mc_st_mc    1.9000e+02 +/-  9.50e-01
  # mean_GausErfExp_tt_mc_tt_mc    1.7374e+02 +/-  6.15e+01
  #         number_st_mc    1.5506e+03 +/-  2.48e+03
  #         number_tt_mc    1.1681e+04 +/-  2.27e+03
  # offset_GausErfExp_tt_mc_tt_mc    1.9999e+02 +/-  1.16e+01
  # sigmaPolyGaus_st_mc_st_mc    1.1456e+01 +/-  1.16e+02
  # sigma_GausErfExp_tt_mc_tt_mc    1.1754e+01 +/-  6.25e+01
  # width_GausErfExp_tt_mc_tt_mc    7.7519e+01 +/-  3.85e+00

  # Puppijet0_msdHi_GausErfExp_tt_mc_tt_mc    3.1671e-01 +/-  3.73e-01
  # Puppijet0_msdHi_PolyGaus_st_mc_st_mc    5.1966e-01 +/-  2.02e+00
  # a1PolyGaus_st_mc_st_mc   -1.9779e-03 +/-  3.01e-03
  # c_GausErfExp_tt_mc_tt_mc   -2.4609e-02 +/-  2.21e-02
  # meanPolyGaus_st_mc_st_mc    1.7878e+02 +/-  9.26e+00
  # mean_GausErfExp_tt_mc_tt_mc    1.7801e+02 +/-  1.36e+00
  #         number_st_mc    1.4015e+03 +/-  2.84e+03
  #         number_tt_mc    1.1830e+04 +/-  2.76e+03
  # offset_GausErfExp_tt_mc_tt_mc    2.0000e+02 +/-  2.88e+00
  # sigmaPolyGaus_st_mc_st_mc    1.8236e+01 +/-  3.63e+01
  # sigma_GausErfExp_tt_mc_tt_mc    1.3885e+01 +/-  2.15e+00
  # width_GausErfExp_tt_mc_tt_mc    7.9302e+01 +/-  1.05e+01

  # converging
  # Puppijet0_msdHi_GausErfExp_tt_mc_tt_mc    1.2856e-01 +/-  5.59e-01
  # Puppijet0_msdHi_PolyGaus_st_mc_st_mc    4.5075e-01 +/-  5.01e-02
  # a1PolyGaus_st_mc_st_mc   -1.9521e-03 +/-  3.21e-03
  # c_GausErfExp_tt_mc_tt_mc   -4.0406e-02 +/-  2.18e-02
  # meanPolyGaus_st_mc_st_mc    1.8042e+02 +/-  2.00e+00
  # mean_GausErfExp_tt_mc_tt_mc    1.7179e+02 +/-  1.82e+01
  #         number_st_mc    7.9511e+03 +/-  3.82e+03
  #         number_tt_mc    5.2802e+03 +/-  3.84e+03
  # offset_GausErfExp_tt_mc_tt_mc    2.0000e+02 +/-  9.70e-01
  # sigmaPolyGaus_st_mc_st_mc    1.4636e+01 +/-  1.15e+01
  # sigma_GausErfExp_tt_mc_tt_mc    9.1787e+00 +/-  1.06e+01
  # width_GausErfExp_tt_mc_tt_mc    5.4555e+01 +/-  1.14e+01

    if "GausExp_tt" in iModel:
        if iLabel == "tt_mc" or iLabel == "tt_data":
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,130.,200.)
        elif "fakeW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,130.,100.,150.)
        elif "realW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,160.,200.)
        lSigma_Gaus        = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        lC_Exp             = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-1.,1.)
        la_Poly           = r.RooRealVar("a_"+lTag,"a_"+lTag,-0.005,0.,0.005)
        la2_Poly           = r.RooRealVar("a2_"+lTag,"a2_"+lTag,-0.005,0.,0.005)
        pGaus              = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pExp               = r.RooExponential("exp_"+lTag,"exp_%s"+lTag,lVar,lC_Exp)
        pPoly              = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(la_Poly,la2_Poly),1)
        lModelPdf          = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pPoly,pExp),r.RooArgList(pVarHigh,pVarHigh1),1)
        if iLabel == 'tt_mc' or iLabel == 'tt_data' or 'realW' in iLabel:
            addConstraint(iWorkspace,lMean_Gaus,175,8,lConstraints)
            addConstraint(iWorkspace,lSigma_Gaus,17,5,lConstraints)
        if 'data' in iModel and 'fakeW' in iLabel:
            addConstraint(iWorkspace,lMean_Gaus,130,10,lConstraints)
            addConstraint(iWorkspace,la_Poly,-0.002,0.001,lConstraints)
            addConstraint(iWorkspace,la2_Poly,-0.005,0.002,lConstraints)

    if iModel == "ErfExp_tt_mc":
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp),r.RooArgList(pVarHigh),1);
        addConstraint(iWorkspace,lMean_Gaus,175,3,lConstraints)
        addConstraint(iWorkspace,lSigma_Gaus,7,2,lConstraints)
        if 'realW' in iLabel:
            addConstraint(iWorkspace,lOffSet_ErfExp,175,10,lConstraints)
            addConstraint(iWorkspace,lWidth_ErfExp,30,5,lConstraints) 

    if "GausGaus_tt" in iModel:
        #narrow component
        lMean_Gaus1     = r.RooRealVar("mean1"+lTag,"mean1"+lTag,175.,140.,200.)
        lSigma_Gaus1    = r.RooRealVar("sigma1"+lTag,"sigma1"+lTag,7.,0.,40.)
        pGaus1          = r.RooGaussian("gaus1"+lTag,"gaus1"+lTag,lVar,lMean_Gaus1,lSigma_Gaus1)
        pVarFrac1       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        pGaus1          = r.RooGaussian("gaus1"+lTag,"gaus1"+lTag,lVar,lMean_Gaus1,lSigma_Gaus1)
        #broad component
        lMean_Gaus2     = r.RooRealVar("mean2"+lTag,"mean2"+lTag,120.,100.,175.) #offset left
        lSigma_Gaus2    = r.RooRealVar("sigma2"+lTag,"sigma2"+lTag,60.,0.,80.) #wider
        pGaus2          = r.RooGaussian("gaus2"+lTag,"gaus2"+lTag,lVar,lMean_Gaus2,lSigma_Gaus2)
        #W
        lMean_Gaus3     = r.RooRealVar("mean3"+lTag,"mean3"+lTag,80.,75.,85.)
        lSigma_Gaus3    = r.RooRealVar("sigma3"+lTag,"sigma3"+lTag,10.,0.,20.)
        pGaus3          = r.RooGaussian("gaus3"+lTag,"gaus3"+lTag,lVar,lMean_Gaus3,lSigma_Gaus3)
        
        #try poly
        a1              = r.RooRealVar("a1"+lTag,"a1"+lTag,0.004,0.,0.007)
        a2              = r.RooRealVar("a2"+lTag,"a2"+lTag,-1e-6,-1e-9,0.)
        a3              = r.RooRealVar("a3"+lTag,"a3"+lTag,1e-7,0.,0.005)
        lExp2             = r.RooPolynomial("lExp2"+lTag,"lExp2"+lTag,lVar,r.RooArgList(a1,a2,a3),1)
        
        #a3              = r.RooRealVar("a3"+lTag,"a3"+lTag,1e-8,-0.005,0.005)
        #lPoly           = r.RooPolynomial("lPoly"+lTag,"lPoly"+lTag,lVar, r.RooArgList(a1,a2),1)
        
        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarFrac1),1);
        #lModelPdf1      = r.RooAddPdf("model_pdf1_%s"%iLabel,"model_pdf1_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarHigh));
        addConstraint(iWorkspace,lMean_Gaus1,175,3,lConstraints)
        addConstraint(iWorkspace,lSigma_Gaus1,7,2,lConstraints)
        lModelPdf       = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,lExp2),r.RooArgList(pVarHigh))

    getattr(iWorkspace,"import")(lModelPdf,r.RooFit.RecycleConflictNodes())
    return iWorkspace.pdf("model_pdf_%s"%iLabel),lConstraints

# return RooExtendPdf and Constraints
def makeModel(iWorkspace,iLabel,iModel):
    print '---- Making model'
    lNumber = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,500,0.,1e7); # check if needed
    lModelPdf,lConstraints = makePdf(iWorkspace,iLabel,iModel)
    lModel = r.RooExtendPdf("model_%s"%iLabel,"model_%s"%iLabel,lModelPdf,lNumber)
    getattr(iWorkspace,"import")(lNumber,r.RooFit.RecycleConflictNodes())
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    iWorkspace.pdf("model_%s"%iLabel).Print()
    return iWorkspace.pdf("model_%s"%iLabel),lConstraints

# fit single rooDataset to model
# return list of constraints
def fitSingleMC(iWorkspace,iLabel,iModel):
    print '---- Fitting for %s'%iLabel
    lVar   = iWorkspace.var(fVar);
    lData  = iWorkspace.data(iLabel+"_D");
    lModel,lConstraints = makeModel(iWorkspace,iLabel,iModel)

    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),                             r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE));
    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"));
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    pRooFitResult.Print()
    getattr(iWorkspace,"import")(pRooFitResult)
    
    draw(lVar,lData,[lModel],pRooFitResult,iLabel+"_only")

    return lConstraints

class TopPeak():
    def __init__(self,options):

        self._lMSD_lo = fXMin
        self._lMSD_hi = fXMax
        self._lMSD_blind_lo = 140;
        self._lMSD_blind_hi = 200;

        self._lMatch  = options.match

        self._lW      = r.RooWorkspace("w","w")

        self._lMSD    = r.RooRealVar(fVar,fVar,self._lMSD_lo,self._lMSD_hi)
        self._lMSD.setRange('Low',self._lMSD_lo,self._lMSD_blind_lo)
        self._lMSD.setRange('Blind',self._lMSD_blind_lo,self._lMSD_blind_hi)
        self._lMSD.setRange('High',self._lMSD_blind_hi,self._lMSD_hi)
        self._lMSD.setBins(fNBins)
        getattr(self._lW,"import")(self._lMSD)

        self._lWeight = r.RooRealVar(fWeight,fWeight,-1e+25,1e+25)
        self._lPt     = r.RooRealVar(fPt,fPt,0,2000)

        self._lPDatas = {}
        self._lHPdfs  = {}
        self._lNPdfs  = {}
        self._lModels = {}
        self._lModels_data= {}
        self._lConstraints = {}

        for iLabel,iFile in fFiles.iteritems():
            pFile = fDir+fTag+"/"+iFile+"_"+fTag+".root"
            if iLabel == 'data': 
                self._lPDatas[iLabel] = self.getDataset(iLabel,pFile,fCut)
            else:
                pCut = fCut
                if 'wlnu' in iLabel: pCut = "(%s&&weight<10)"%fCut # get rid of high-weight event for WLNu
                self._lPDatas[iLabel] = self.getDataset(iLabel,pFile,pCut,True)
                # check single fits first
                if iLabel != "mc":
                    self._lConstraints[iLabel] = fitSingleMC(self._lW,iLabel,fPDFs[iLabel])

            # check why are these needed
            self._lHPdfs[iLabel] = self.histPdf(iLabel)
            self._lNPdfs[iLabel] = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,self._lW.data(iLabel+"_D").sumEntries()) # *sf? ttsf, what is that?
            self._lNPdfs[iLabel].Print()
            print self._lW.data(iLabel+"_D").sumEntries()
        
    # get dataset for each process from tree
    def getDataset(self,iLabel,iFile,iCut="(1==1)",iWeight=False,iTree="otree2"):
        print 'preparing dataset for ',iFile,' with cut ',iCut
        lFile   = r.TFile(iFile)
        lTree   = lFile.Get(iTree)
        #lData   = r.RooDataSet(iLabel+"_D",iLabel+"_D",lTree,r.RooAet(self._lMSD,self._lPt,self._lWeight),iCut)
        lData = r.RooDataSet(iLabel+"_D",iLabel+"_D",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight))
        lCut = r.TTreeFormula("lCut",iCut,lTree)  
        lTree.SetNotify(lCut)       
        for i0 in range(lTree.GetEntriesFast()):
        #for i0 in range(1000): 
            lTree.LoadTree(i0) 
            lSel = False 
            for i1 in range(lCut.GetNdata()):  
                if (lCut.EvalInstance(i1)):  
                    lSel=True; break;
            if not lSel: continue     
            lTree.GetEntry(i0) 
            if iWeight:
                lWeight = getattr(lTree,fWeight)*fLumi
            else:
                lWeight = 1
            lMass = getattr(lTree,fVar)
            lMatched = 0
            jmatched = getattr(lTree,"Puppijet0_vMatching");
            jhadronic = getattr(lTree,"Puppijet0_isHadronicV");
            if jhadronic == 1.0 and jmatched < 1. and jmatched > 0.:
                lMatched = 1;
            if 'realW' in iLabel and lMatched == 0: continue
            if 'fakeW' in iLabel and lMatched == 1: continue
            if lMass < self._lMSD_hi and lMass > self._lMSD_lo:
                self._lMSD.setVal(lMass)
                lData.add(r.RooArgSet(self._lMSD), lWeight)
        getattr(self._lW,'import')(lData,r.RooFit.RecycleConflictNodes())
        lData.Print()
        lFile.Close()
        return self._lW.data(iLabel+"_D")
                     
    # shape and normalization parameters fixed
    def getModel(self, iLabel):
        print 'getting model for ', iLabel
        lModel    = self._lW.pdf("model_"+iLabel)
        lModel.Print()
        self._lW.data(iLabel+"_D").Print()
        # lPars     = lModel.getParameters(self._lW.data(iLabel+"_D"))
        # pParsIter = lPars.createIterator()
        # pParsIter.Reset()
        # pPar = pParsIter.Next()
        # while (pPar):
        #     pPar.setConstant(r.kTRUE)
        #     pPar = pParsIter.Next()
        return lModel

    # return roohistpdf
    def histPdf(self, iLabel):
        lDataSet = self._lW.data(iLabel+"_D")
        lVar = self._lW.var(fVar)
        lReduced  = lDataSet.reduce(r.RooArgSet(lVar))
        lDataHist = lReduced.binnedClone(lReduced.GetName()+"_binnedClone",lReduced.GetName()+"_binnedClone");
        lHPdf     = r.RooHistPdf(lReduced.GetName()+"P",lReduced.GetName()+"P",r.RooArgSet(lVar),lDataHist,0)
        print lReduced.sumEntries()
        getattr(self._lW,"import")(lHPdf,r.RooFit.RecycleConflictNodes())
        return lHPdf

    # fit mc and data
    def fit(self):

        lVar   = self._lW.var(fVar);
        lSData = self._lPDatas['data']
        lSMc   = self._lPDatas['mc']

        # constraints
        pPdfConstraints_Mc   = r.RooArgSet("pPdfConstraints_Mc")
        pPdfConstraints_Data = r.RooArgSet("pPdfConstraints_Data")
        pSamples = ['st_mc','st_data','tt_mc','tt_data']
        if self._lMatch:
            pSamples = ['st_mc','st_data','tt_fakeW_mc','tt_fakeW_data','tt_realW_mc','tt_realW_data']
        for iLabel in pSamples:
            self._lModels[iLabel] = self.getModel(iLabel)
            self._lModels[iLabel].Print()
            for i0 in range(len(self._lConstraints[iLabel])):
                if 'mc' in iLabel:
                    pPdfConstraints_Mc.add(self._lW.pdf(self._lConstraints[iLabel][i0]))
                if 'data' in iLabel:
                    pPdfConstraints_Data.add(self._lW.pdf(self._lConstraints[iLabel][i0]))
        pPdfConstraints_Mc.Print()
        pPdfConstraints_Data.Print()

        # models
        if self._lMatch:
            self._lModels['Mc'] = r.RooAddPdf("model_mc","model_mc",r.RooArgList(self._lModels['st_mc'],self._lModels['tt_fakeW_mc'],self._lModels['tt_realW_mc']))
            self._lModels['Data'] = r.RooAddPdf("model_data","model_data",r.RooArgList(self._lModels['st_data'],self._lModels['tt_fakeW_data'],self._lModels['tt_realW_data']))
        else:
            self._lModels['Mc'] = r.RooAddPdf("model_mc","model_mc",r.RooArgList(self._lModels['st_mc'],self._lModels['tt_mc']))
            self._lModels['Data'] = r.RooAddPdf("model_data","model_data",r.RooArgList(self._lModels['st_data'],self._lModels['tt_data']))

        # fit to mc
        pRooFitResult_Mc   = self._lModels['Mc'].fitTo(lSMc,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Mc),r.RooFit.Verbose(r.kFALSE))
        pRooFitResult_Mc   = self._lModels['Mc'].fitTo(lSMc,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Mc))
        pRooFitResult_Mc.Print()
        lTagMc = 'mc_only'
        if self._lMatch: lTagMc += '_ttmatched'
        draw(lVar,lSMc,[self._lModels['Mc']],pRooFitResult_Mc,lTagMc)
 
        # print mc parameters
        x1 = self._lModels['Mc'].getVariables()
        print x1.Print("v")

        # fit to data
        pRooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Data),r.RooFit.Verbose(r.kFALSE))
        pRooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.ExternalConstraints(pPdfConstraints_Data))
        pRooFitResult_Data.Print()
        lTagData = 'data_only'
        if self._lMatch: lTagData+= '_ttmatched'
        draw(lVar,lSData,[self._lModels['Data']],pRooFitResult_Data,lTagData) 

        # print data parameters
        x2 = self._lModels['Data'].getVariables()
        print x2.Print("v")

        params = {}
        if self._lMatch:
            lTagMc = fPDFs['tt_realW_mc']+'_tt_realW_mc'
            lTagData = fPDFs['tt_realW_data']+'_tt_realW_data'
        else:
            lTagMc = fPDFs['tt_mc']+'_tt_mc'
            lTagData = fPDFs['tt_data']+'_tt_data'
        params["mc_mean_val"] = x1.find("mean_%s"%lTagMc).getVal()
        params["mc_mean_err"] = x1.find("mean_%s"%lTagMc).getError()
        params["mc_sigma_val"]   = x1.find("sigma_%s"%lTagMc).getVal() 
        params["mc_sigma_err"]   = x1.find("sigma_%s"%lTagMc).getError()

        params["data_mean_val"] = x2.find("mean_%s"%lTagData).getVal()
        params["data_mean_err"] = x2.find("mean_%s"%lTagData).getError()
        params["data_sigma_val"]   = x2.find("sigma_%s"%lTagData).getVal()
        params["data_sigma_err"]   = x2.find("sigma_%s"%lTagData).getError()
        print params

        # draw data and mc
        lTagDataMc = "data_mc"
        if self._lMatch: lTagDataMc += "_ttmatched"
        drawDataMc(lVar,lSData,[self._lModels['Data']],lSMc,[self._lModels['Mc']],pRooFitResult_Mc,pRooFitResult_Data,params,lTagDataMc)
        # write workspace
        self._lW.writeToFile(fOutput)


if __name__ == "__main__":
    options = parser()
    print options
    global fTag,fXMin,fXMax,fNBins
    fTag   = options.tag
    fXMin  = options.xmin
    fXMax  = options.xmax
    fNBins = int( (fXMax - fXMin) / fBinWidth)
    # get roodataset and make individual fits
    lTop = TopPeak(options);
    # make combined fit
    lTop.fit();

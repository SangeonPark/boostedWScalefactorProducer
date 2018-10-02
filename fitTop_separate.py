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
          #'vv':   'VV', # ignore diboson for now
          'st_mc':   'ST',
          'st_data':   'ST',
          'tt_mc':   'TT',
          'tt_data':   'TT',
          #'tt_fakeW': 'TT',
          #'tt_realW': 'TT',
          'mc':   'PseudoData'}

fPDFs = {
         #'wlnu_mc': 'Exp_mc',
         #'wlnu_data': 'Exp_data',
         #'vv':   '',
         #'st_mc':    'GausErfExp_ttbar',
         'st_mc':   'PolyGaus_st_mc',
         'st_data':   'PolyGaus_st_data',
         #'tt':   'GausErfExp_ttbar',
         'tt_mc':    'GausGaus_ttbar_mc',
         'tt_data':    'GausGaus_ttbar_data',
         #'tt_fakeW': 'GausGaus_ttbar', #what is the real shape??????
         #'tt_realW': 'GausGaus_ttbar' #what is the real shape??????
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

def draw(iVar,iData,iFuncs,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    drawFrame(lFrame,iData,iFuncs)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %f GeV"%fBinWidth);
    lFrame.Draw()
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
    print iData.GetName()+"_binnedClone"
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
def addConstraint(iVar,iMean,iSigma,iList):
    lMean = r.RooRealVar("%s_mean"%iVar.GetName(),"%s_mean"%iVar.GetName(),iMean);
    lSigma = r.RooRealVar("%s_sigma"%iVar.GetName(),"%s_sigma"%iVar.GetName(),iSigma);
    lConstraint = r.RooGaussian("constraint_pdf_%s"%iVar.GetName(),"constraint_pdf_%s"%iVar.GetName(),iVar,lMean,lSigma)
    iList.append(lConstraint.GetName())
    return lConstraint

# make Pdf from model, probably should put parameters in dictionary
def makePdf(iWorkspace,iLabel,iModel,iMc=False):
    lVar = iWorkspace.var(fVar);        
    print 'making pdf' 
    lModelPdf = None
    lTag = "%s_%s"%(iModel,iLabel)
    pVarHigh       = r.RooRealVar("%sHi_%s"%(lVar.GetName(),lTag),"%sHi_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    pVarHigh1      = r.RooRealVar("%sHi1_%s"%(lVar.GetName(),lTag),"%sHi1_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    
    if iModel == "Exp_mc":
        lC_Exp_mc         = r.RooRealVar("c_mc"+lTag,"c_mc"+lTag,-0.03,-2.,0.05)
        lModelPdf      = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_mc);
        
    if iModel == "Exp_data":
        lC_Exp_data         = r.RooRealVar("c_data"+lTag,"c_data"+lTag,-0.03,-2.,0.05)
        lModelPdf      = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_data);
        
    if iModel == "PolyGaus_st_mc":
        lMean_Gaus_mc      = r.RooRealVar("mean_mc"+lTag,"mean_mc"+lTag,175.,120.,190.)
        lSigma_Gaus_mc     = r.RooRealVar("sigma_mc"+lTag,"sigma_mc"+lTag,15.,0.,100.)
        pGaus_mc           = r.RooGaussian("gaus_mc"+lTag,"gaus_mc%s"+lTag,lVar,lMean_Gaus_mc,lSigma_Gaus_mc)
        lMean2_Gaus_mc      = r.RooRealVar("mean2_mc"+lTag,"mean2_mc"+lTag,150.,140.,175.)
        lSigma2_Gaus_mc     = r.RooRealVar("sigma2_mc"+lTag,"sigma2_mc"+lTag,80.,40.,100.)
        pGaus2_mc           = r.RooGaussian("gaus2_mc"+lTag,"gaus2_mc%s"+lTag,lVar,lMean2_Gaus_mc,lSigma2_Gaus_mc)
        la1_mc             = r.RooRealVar("a1_mc"+lTag,"a1_mc"+lTag,-0.003,-0.005,0.)
        #la1_mc             = r.RooRealVar("a1_mc"+lTag,"a1_mc"+lTag,-0.04,-.1,0.)
        #la2_mc             = r.RooRealVar("a2_mc"+lTag,"a2_mc"+lTag,0.1,0.,1.)
        #pPoly_mc           = r.RooExponential("poly_mc"+lTag,"poly_mc%s"+lTag,lVar,la1_mc)
        pPoly_mc           = r.RooPolynomial("poly_mc"+lTag,"poly_mc%s"+lTag,lVar,r.RooArgList(la1_mc),1)
        lModelPdf          = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus_mc,pPoly_mc),r.RooArgList(pVarHigh))

    if iModel == "PolyGaus_st_data":
        lMean_Gaus_data      = r.RooRealVar("mean_data"+lTag,"mean_data"+lTag,175.,120.,190.)
        lSigma_Gaus_data     = r.RooRealVar("sigma_data"+lTag,"sigma_data"+lTag,15.,0.,100.)
        pGaus_data           = r.RooGaussian("gaus_data"+lTag,"gaus_data%s"+lTag,lVar,lMean_Gaus_data,lSigma_Gaus_data)
        #lMean2_Gaus_data      = r.RooRealVar("mean2_data"+lTag,"mean2_data"+lTag,150.,140.,175.)
        #lSigma2_Gaus_data     = r.RooRealVar("sigma2_data"+lTag,"sigma2_data"+lTag,80.,40.,100.)
        #pGaus2_data           = r.RooGaussian("gaus2_data"+lTag,"gaus2_data%s"+lTag,lVar,lMean2_Gaus_data,lSigma2_Gaus_data)
        la1_data             = r.RooRealVar("a1_data"+lTag,"a1_data"+lTag,-0.003,-0.005,0.)
        #la1_data             = r.RooRealVar("a1_data"+lTag,"a1_data"+lTag,-0.04,-.1,0.)
        #la2_data             = r.RooRealVar("a2_data"+lTag,"a2_data"+lTag,0.1,0.,1.)
        pPoly_data             = r.RooPolynomial("poly_data"+lTag,"poly_data%s"+lTag,lVar,r.RooArgList(la1_data),1)
        #pPoly_data           = r.RooExponential("poly_data"+lTag,"poly_data%s"+lTag,lVar,la1_data)
        lModelPdf            = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus_data,pPoly_data),r.RooArgList(pVarHigh))

    if iModel == "ErfExpGaus_st_mc":
        lC_ErfExp_mc      = r.RooRealVar("c_mc"+lTag,"c_mc"+lTag,-0.05,-0.2,0.0)
        lOffSet_ErfExp_mc  = r.RooRealVar("offset_mc"+lTag,"offset_mc"+lTag,175.,140.,300.)
        lWidth_ErfExp_mc   = r.RooRealVar("width_mc"+lTag,"width_mc"+lTag,30.,10.,300.)
        lMean_Gaus_mc      = r.RooRealVar("mean_mc"+lTag,"mean_mc"+lTag,175.,140.,200.)
        lSigma_Gaus_mc     = r.RooRealVar("sigma_mc"+lTag,"sigma_mc"+lTag,7.,0.,40.)
        pGaus_mc           = r.RooGaussian("gaus_mc"+lTag,"gaus_mc%s"+lTag,lVar,lMean_Gaus_mc,lSigma_Gaus_mc)
        #add W 
        #lMean_Gaus_mc2    = r.RooRealVar("mean_mc2"+lTag,"mean_mc2"+lTag,85.,75.,95.)
        #lSigma_Gaus_mc2   = r.RooRealVar("sigma_mc2"+lTag,"sigma_mc2"+lTag,10.,2.,20.)
        #pGaus_mc2         = r.RooGaussian("gaus_mc2"+lTag,"gaus_mc2%s"+lTag,lVar,lMean_Gaus_mc2,lSigma_Gaus_mc2)
        
        pErfExp_mc         = r.RooErfExpPdf("erfExp_mc"+lTag,"erfExp_mc%s"+lTag,lVar,lC_ErfExp_mc,lOffSet_ErfExp_mc,lWidth_ErfExp_mc);
        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp_mc,pGaus_mc,pGaus_mc2),r.RooArgList(pVarHigh,pVarHigh1));
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp_mc,pGaus_mc),r.RooArgList(pVarHigh));

    if iModel == "ErfExpGaus_st_data":
        lC_ErfExp_data      = r.RooRealVar("c_data"+lTag,"c_data"+lTag,-0.05,-0.2,0.0)
        lOffSet_ErfExp_data = r.RooRealVar("offset_data"+lTag,"offset_data"+lTag,175.,140.,300.)
        lWidth_ErfExp_data  = r.RooRealVar("width_data"+lTag,"width_data"+lTag,30.,10.,300.)
        lMean_Gaus_data     = r.RooRealVar("mean_data"+lTag,"mean_data"+lTag,175.,140.,200.)
        lSigma_Gaus_data    = r.RooRealVar("sigma_data"+lTag,"sigma_data"+lTag,7.,0.,40.)
        pGaus_data          = r.RooGaussian("gaus_data"+lTag,"gaus_data%s"+lTag,lVar,lMean_Gaus_data,lSigma_Gaus_data)
        #add W 
        #lMean_Gaus_data2    = r.RooRealVar("mean_data2"+lTag,"mean_data2"+lTag,85.,75.,95.)
        #lSigma_Gaus_data2   = r.RooRealVar("sigma_data2"+lTag,"sigma_data2"+lTag,10.,2.,20.)
        #pGaus_data2         = r.RooGaussian("gaus_data2"+lTag,"gaus_data2%s"+lTag,lVar,lMean_Gaus_data2,lSigma_Gaus_data2)

        pErfExp_data        = r.RooErfExpPdf("erfExp_data"+lTag,"erfExp_data%s"+lTag,lVar,lC_ErfExp_data,lOffSet_ErfExp_data,lWidth_ErfExp_data);
        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp_data,pGaus_data,pGaus_data2),r.RooArgList(pVarHigh,pVarHigh1));
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp_data,pGaus_data),r.RooArgList(pVarHigh));

    if iModel == "GausErfExp_ttbar":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,0.,200.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,30.,0.,200.)
        lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,130.,200.)
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        pErfExp        = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_%s"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pVarFrac       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        lC_ErfExp.setConstant(r.kTRUE);
        lOffSet_ErfExp.setConstant(r.kTRUE);
        lWidth_ErfExp.setConstant(r.kTRUE);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp),r.RooArgList(pVarHigh),1);
        
    if iModel == "ErfExp_st_mc":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf      = pErfExp
        
    if iModel == "ErfExp_st_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf      = pErfExp
        
    if iModel == "GausGaus_ttbar_mc":
        #narrow component
        lMean_Gaus1_mc     = r.RooRealVar("mean1_mc"+lTag,"mean1_mc"+lTag,175.,140.,200.)
        lSigma_Gaus1_mc    = r.RooRealVar("sigma1_mc"+lTag,"sigma1_mc"+lTag,7.,0.,40.)
        pGaus1_mc          = r.RooGaussian("gaus1_mc"+lTag,"gaus1_mc%s"+lTag,lVar,lMean_Gaus1_mc,lSigma_Gaus1_mc)
        pVarFrac1_mc       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        pGaus1_mc          = r.RooGaussian("gaus1_mc"+lTag,"gaus1_mc%s"+lTag,lVar,lMean_Gaus1_mc,lSigma_Gaus1_mc)
        #broad component
        lMean_Gaus2_mc     = r.RooRealVar("mean2_mc"+lTag,"mean2_mc"+lTag,120.,100.,175.) #offset left
        lSigma_Gaus2_mc    = r.RooRealVar("sigma2_mc"+lTag,"sigma2_mc"+lTag,60.,0.,80.) #wider
        pGaus2_mc          = r.RooGaussian("gaus2_mc"+lTag,"gaus2_mc%s"+lTag,lVar,lMean_Gaus2_mc,lSigma_Gaus2_mc)
        #W
        lMean_Gaus3_mc     = r.RooRealVar("mean3_mc"+lTag,"mean3_mc"+lTag,80.,75.,85.)
        lSigma_Gaus3_mc    = r.RooRealVar("sigma3_mc"+lTag,"sigma3_mc"+lTag,10.,0.,20.)
        pGaus3_mc          = r.RooGaussian("gaus3_mc"+lTag,"gaus3_mc%s"+lTag,lVar,lMean_Gaus3_mc,lSigma_Gaus3_mc)
        
        #try poly
        a1_mc              = r.RooRealVar("a1_mc"+lTag,"a1_mc"+lTag,0.004,0.,0.007)
        a2_mc              = r.RooRealVar("a2_mc"+lTag,"a2_mc"+lTag,-1e-6,-1e-9,0.)
        a3_mc              = r.RooRealVar("a3_mc"+lTag,"a3_mc"+lTag,1e-7,0.,0.005)
        lExp2_mc             = r.RooPolynomial("lExp2_mc"+lTag,"lExp2_mc"+lTag,lVar,r.RooArgList(a1_mc,a2_mc,a3_mc),1)
        
        #a3_mc              = r.RooRealVar("a3_mc"+lTag,"a3_mc"+lTag,1e-8,-0.005,0.005)
        #lPoly_mc           = r.RooPolynomial("lPoly_mc"+lTag,"lPoly_mc"+lTag,lVar, r.RooArgList(a1_mc,a2_mc),1)
        
        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarFrac1),1);
        #lModelPdf1      = r.RooAddPdf("model_pdf1_%s"%iLabel,"model_pdf1_%s"%iLabel,r.RooArgList(pGaus1_mc,pGaus2_mc),r.RooArgList(pVarHigh));
        lModelPdf       = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1_mc,lExp2_mc),r.RooArgList(pVarHigh))

    if iModel == "GausGaus_ttbar_data":
        #narrow component
        lMean_Gaus1_data     = r.RooRealVar("mean1_data"+lTag,"mean1_data"+lTag,175.,140.,200.)
        lSigma_Gaus1_data    = r.RooRealVar("sigma1_data"+lTag,"sigma1_data"+lTag,7.,0.,30.)
        pGaus1_data          = r.RooGaussian("gaus1_data"+lTag,"gaus1_data%s"+lTag,lVar,lMean_Gaus1_data,lSigma_Gaus1_data)
        pVarFrac1_data       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        pGaus1_data          = r.RooGaussian("gaus1_data"+lTag,"gaus1_data%s"+lTag,lVar,lMean_Gaus1_data,lSigma_Gaus1_data)
        #broad component
        lMean_Gaus2_data     = r.RooRealVar("mean2_data"+lTag,"mean2_data"+lTag,120.,100.,175.) #offset left
        lSigma_Gaus2_data    = r.RooRealVar("sigma2_data"+lTag,"sigma2_data"+lTag,60.,0.,80.) #wider
        pGaus2_data          = r.RooGaussian("gaus2_data"+lTag,"gaus2_data%s"+lTag,lVar,lMean_Gaus2_data,lSigma_Gaus2_data)
        #W
        lMean_Gaus3_data     = r.RooRealVar("mean3_data"+lTag,"mean3_data"+lTag,80.,75.,85.)
        lSigma_Gaus3_data    = r.RooRealVar("sigma3_data"+lTag,"sigma3_data"+lTag,10.,0.,20.)
        pGaus3_data          = r.RooGaussian("gaus3_data"+lTag,"gaus3_data%s"+lTag,lVar,lMean_Gaus3_data,lSigma_Gaus3_data)
        #try poly
        a1_data              = r.RooRealVar("a1_data"+lTag,"a1_data"+lTag,0.004,0.,0.007)
        a2_data              = r.RooRealVar("a2_data"+lTag,"a2_data"+lTag,-1e-6,-1e-9,0.)  
        a3_data              = r.RooRealVar("a3_data"+lTag,"a3_data"+lTag,1e-7,0.,0.005)
        lPoly_data           = r.RooPolynomial("lPoly_data"+lTag,"lPoly_data"+lTag,lVar,r.RooArgList(a1_data,a2_data,a3_data),1)

        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarFrac1),1);
        #lModelPdf1      = r.RooAddPdf("model_pdf1_%s"%iLabel,"model_pdf1_%s"%iLabel,r.RooArgList(pGaus1_data,pGaus2_data),r.RooArgList(pVarHigh));
        lModelPdf       = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1_data,lPoly_data),r.RooArgList(pVarHigh))

    getattr(iWorkspace,"import")(lModelPdf,r.RooFit.RecycleConflictNodes())
    return iWorkspace.pdf("model_pdf_%s"%iLabel)

# return RooExtendPdf
def makeModel(iWorkspace,iLabel,iModel):
    print 'making model'
    lNumber = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,500,0.,1e7); # check if needed
    lModelPdf = makePdf(iWorkspace,iLabel,iModel)
    lModelPdf.Print()
    lModel = r.RooExtendPdf("model_%s"%iLabel,"model_%s"%iLabel,lModelPdf,lNumber)
    getattr(iWorkspace,"import")(lNumber,r.RooFit.RecycleConflictNodes())
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    lModel.Print();
    iWorkspace.pdf("model_%s"%iLabel).Print()
    return iWorkspace.pdf("model_%s"%iLabel)

# fit single rooDataset to model
def fitSingleMC(iWorkspace,iLabel,iModel):
    print '---fitting for %s'%iLabel
    lVar   = iWorkspace.var(fVar);
    lData  = iWorkspace.data(iLabel+"_D");
    lModel = makeModel(iWorkspace,iLabel,iModel)

    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),                             r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE));
    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"));
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    pRooFitResult.Print()
    getattr(iWorkspace,"import")(pRooFitResult)

    draw(lVar,lData,[lModel],iLabel+"_only")

class TopPeak():
    def __init__(self,options):

        self._lMSD_lo = fXMin
        self._lMSD_hi = fXMax
        self._lMSD_blind_lo = 140;
        self._lMSD_blind_hi = 200;

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
                    fitSingleMC(self._lW,iLabel,fPDFs[iLabel])

            # check why are these needed
            self._lHPdfs[iLabel] = self.histPdf(iLabel)
            print self._lPDatas[iLabel]
            print self._lW.data(iLabel+"_D")
            self._lNPdfs[iLabel] = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,self._lW.data(iLabel+"_D").sumEntries()) # *sf? ttsf, what is that?
            print 'number in pdf ',self._lNPdfs[iLabel]
            self._lNPdfs[iLabel].Print()
            print self._lW.data(iLabel+"_D").sumEntries()

    # get dataset for each process from tree
    def getDataset(self,iLabel,iFile,iCut="(1==1)",iWeight=False,iTree="otree2"):
        print 'preparing dataset for ',iFile,' with cut ',iCut
        lFile   = r.TFile(iFile)
        lTree   = lFile.Get(iTree)
        print lTree
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
            #TODO add separate samples for tt_fakeW and tt_realW. if tt_fakeW, reject if LeadingJet_MatchedHadW = True and vice versa
            #if iLabel == "tt_realW" and getattr(lTree,"LeadingJet_MatchedHadW") == 0.0: continue
            #if iLabel == "tt_fakeW" and getattr(lTree,"LeadingJet_MatchedHadW") == 1.0: continue 
            if iWeight:
                lWeight = getattr(lTree,fWeight)*fLumi
            else:
                lWeight = 1
            lMass   = getattr(lTree,fVar) 
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
        # create extended pdf
        #lNorm = r.RooRealVar("norm_"+iLabel,"norm_"+iLabel,500.,0.,1e5);
        #lModelExtend = r.RooExtendPdf("model_"+iLabel,"model_"+iLabel,lModel,lNorm)
        #getattr(self._lW,'import')(lModelExtend,r.RooFit.RecycleConflictNodes())
        #return lModelExtend
        return lModel

    # retunr roohistpdf
    def histPdf(self, iLabel):
        lDataSet = self._lW.data(iLabel+"_D")
        lVar = self._lW.var(fVar)
        lReduced = lDataSet.reduce(r.RooArgSet(lVar))
        lDataHist = lReduced.binnedClone(lReduced.GetName()+"_binnedClone",lReduced.GetName()+"_binnedClone");
        lDataHist.Print()
        lHPdf     = r.RooHistPdf(lReduced.GetName()+"P",lReduced.GetName()+"P",r.RooArgSet(lVar),lDataHist,0)
        lHPdf.Print()
        print lReduced.sumEntries()
        getattr(self._lW,"import")(lHPdf,r.RooFit.RecycleConflictNodes())
        return lHPdf

    # fit mc and data
    def fit(self):

        lVar   = self._lW.var(fVar);
        lSData = self._lPDatas['data']
        lSMc   = self._lPDatas['mc']
        print 'DATA', lSData.Print()

        print 'MODELS'
        #print self._lModels
        #for iLabel in ['tt','st','wlnu']:
        for iLabel in ['st_mc','st_data','tt_mc','tt_data']:
            self._lModels[iLabel] = self.getModel(iLabel)
            self._lModels[iLabel].Print()
        '''
        # self._lModels['Mc']  = r.RooAddPdf("model_mc","model_mc",
        #                                    r.RooArgList(self._lHPdfs['tt'],self._lHPdfs['st'],self._lHPdfs['wlnu']),
        #                                    r.RooArgList(self._lNPdfs['st'],self._lNPdfs['wlnu']))
        # so apparatenly there should be 2 type of "ttbar" models: model_bkg_data,model_ttbar_data, and also: model_bkg_TotalMC, model_ttbar_TotalMC and I do not understnd the differences...
        self._lModels['Mc'] = r.RooAddPdf("model_mc","model_mc",r.RooArgList(self._lModels['st_mc'], self._lModels['tt_mc']))
        self._lModels['Data'] = r.RooAddPdf("model_data","model_data",r.RooArgList(self._lModels['st_data'],self._lModels['tt_data']))

        #getattr(self._lW,'import')(self._lModels['Mc'],r.RooFit.RecycleConflictNodes())
        
        #print self._lModels['tt']
        # fit to mc
        pRooFitResult_Mc   = self._lModels['Mc'].fitTo(lSMc,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
        pRooFitResult_Mc.Print()
        draw(lVar,lSMc,[self._lModels['Mc']],'mc_only')
 
        x1 = self._lModels['Mc'].getVariables()
        print x1.Print("v")

        #fit to data
        pRooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Strategy(1),r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE)) #r.ExternalConstraints(lDataConstraints))
        pRooFitResult_Data.Print()
        draw(lVar,lSData,[self._lModels['Data']],'data_only') # make this another frame? with chi2? and include pull?

        x2 = self._lModels['Data'].getVariables()
        print x2.Print("v")

        params = {}
        params["mc_mean_val"] = x1.find("mean1_mcGausGaus_ttbar_mc_tt_mc").getVal()
        params["mc_mean_err"] = x1.find("mean1_mcGausGaus_ttbar_mc_tt_mc").getError()
        params["mc_sigma_val"]   = x1.find("sigma1_mcGausGaus_ttbar_mc_tt_mc").getVal() 
        params["mc_sigma_err"]   = x1.find("sigma1_mcGausGaus_ttbar_mc_tt_mc").getError()

        params["data_mean_val"] = x2.find("mean1_dataGausGaus_ttbar_data_tt_data").getVal()
        params["data_mean_err"] = x2.find("mean1_dataGausGaus_ttbar_data_tt_data").getError()
        params["data_sigma_val"]   = x2.find("sigma1_dataGausGaus_ttbar_data_tt_data").getVal()
        params["data_sigma_err"]   = x2.find("sigma1_dataGausGaus_ttbar_data_tt_data").getError()

        print params

        # draw data and mc
        drawDataMc(lVar,lSData,[self._lModels['Data']],lSMc,[self._lModels['Mc']],pRooFitResult_Mc,pRooFitResult_Data,params,"data_mc")
        '''
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
    # get roodataset
    lTop = TopPeak(options);
    # make combined fit
    lTop.fit();

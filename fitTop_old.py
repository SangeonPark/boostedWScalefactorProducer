#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser

import tdrstyle
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
          'wlnu': 'WJets',
          #'vv':   'VV', # ignore diboson for now
          'st':   'ST',
          'tt':   'TT',
          'mc':   'PseudoData'}

fPDFs = {
         'wlnu': 'Exp',
         #'vv':   '',
         'st':   'ErfExpGaus_sp',
         #'tt':   'GausErfExp_ttbar',
         'tt':    'GausGaus_ttbar',
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

def drawDataMc(iVar,iData,iFuncsData,iMc,iFuncsMc,iLabel="A"):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    drawFrame(lFrame,iData,iFuncsData,True)
    drawFrame(lFrame,iMc,iFuncsMc)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %f GeV"%fBinWidth);
    lFrame.Draw()
    lCan.Modified()
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")

def getChi2NDOF(iData,iModel,iRooFitResult,iVar):
    lDataHist = iData.binnedClone(iData.GetName()+"_binnedClone",iData.GetName()+"_binnedClone");
    pChi2     = iModel.createChi2(lDataHist,r.RooFit.Extended(r.kTRUE),r.RooFit.DataError(r.RooAbsData.Poisson))
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

# make Pdf from model, probably should put parameters in dictionary
def makePdf(iWorkspace,iLabel,iModel,iMc=False):
    lVar = iWorkspace.var(fVar);        
    pVarHigh       = r.RooRealVar("%sHi_%s"%(lVar.GetName(),iLabel),"%sHi_%s"%(lVar.GetName(),iLabel),0.5,0.,1.);

    lModelPdf = None
    lTag = "%s_%s"%(iModel,iLabel)

    if iModel == "Exp":
        lC_Exp         = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.03,-2.,0.05)
        lModelPdf      = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp);

    if iModel == "ErfExpGaus_sp":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.04,-0.2,0.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,140.,200.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,30.,10.,300.)
        lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,140.,200.)
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_%s"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus),r.RooArgList(pVarHigh));
        
    # isn't this the exact same pdf as above?
    if iModel == "GausErfExp_ttbar":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,0.,200.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,30.,0.,200.)
        lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,175.,140.,200.)
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        pErfExp        = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_%s"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pVarFrac       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        lC_ErfExp.setConstant(r.kTRUE);
        lOffSet_ErfExp.setConstant(r.kTRUE);
        lWidth_ErfExp.setConstant(r.kTRUE);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp),r.RooArgList(pVarFrac),1);

    if iModel == "ErfExp_ttbar":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,0.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,175.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf      = pErfExp
        # there should be constraints here?
#         lGaus1         = r.RooGaussian.addConstraint(lC_ErfExp,lC_ErfExp.getVal(),-0.02,constraint)
#         getattr(iWorkspace,"import")(lGaus1,r.RooFit.RecycleConflictNodes())

    if iModel == "GausGaus_ttbar":
        #narrow component
        lMean_Gaus1     = r.RooRealVar("mean1_"+lTag,"mean1_"+lTag,175.,140.,200.)
        lSigma_Gaus1    = r.RooRealVar("sigma1_"+lTag,"sigma1_"+lTag,7.,0.,20.)
        pGaus1          = r.RooGaussian("gaus1_"+lTag,"gaus1_%s"+lTag,lVar,lMean_Gaus1,lSigma_Gaus1)
        pVarFrac1       = r.RooRealVar("%sFrac_%s"%(lVar.GetName(),iLabel),"%sFrac_%s"%(lVar.GetName(),iLabel),1.0); # is this needed?
        pGaus1          = r.RooGaussian("gaus1_"+lTag,"gaus1_%s"+lTag,lVar,lMean_Gaus1,lSigma_Gaus1)
        #broad component
        lMean_Gaus2     = r.RooRealVar("mean2_"+lTag,"mean2_"+lTag,120.,100.,175.) #offset left
        lSigma_Gaus2    = r.RooRealVar("sigma2_"+lTag,"sigma2_"+lTag,60.,0.,80.) #wider
        pGaus2          = r.RooGaussian("gaus2_"+lTag,"gaus2_%s"+lTag,lVar,lMean_Gaus2,lSigma_Gaus2)
        
        #lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarFrac1),1);
        lModelPdf      = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus1,pGaus2),r.RooArgList(pVarHigh));

    getattr(iWorkspace,"import")(lModelPdf,r.RooFit.RecycleConflictNodes())
    return iWorkspace.pdf("model_pdf_%s"%iLabel)

# return RooExtendPdf
def makeModel(iWorkspace,iLabel,iModel):
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
        self._lConstraints = {}

        for iLabel,iFile in fFiles.iteritems():
            pFile = fDir+fTag+"/"+iFile+"_"+fTag+".root"
            if 'data' in iLabel: 
                self._lPDatas[iLabel] = self.getDataset(iLabel,pFile,fCut)
            else:
                pCut = fCut
                if 'wlnu' in iLabel: pCut = "(%s&&weight<10)"%fCut # get rid of high-weight event for WLNu
                pCut = "%f*weight*(%s)"%(fLumi,pCut)
                self._lPDatas[iLabel] = self.getDataset(iLabel,pFile,pCut)
                
                # check single fits first
                if 'mc' not in iLabel:
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
    def getDataset(self,iLabel,iFile,iCut="(1==1)",iTree="otree2"):
        print 'preparing dataset for ',iFile
        lFile   = r.TFile(iFile)
        lTree   = lFile.Get(iTree)
        print lTree
        lData   = r.RooDataSet(iLabel+"_D",iLabel+"_D",lTree,r.RooArgSet(self._lMSD,self._lPt,self._lWeight),iCut)
        getattr(self._lW,'import')(lData,r.RooFit.RecycleConflictNodes())
        lData.Print()
        lFile.Close()
        return self._lW.data(iLabel+"_D")
                     
    # shape and normalization parameters fixed
    def getModel(self, iLabel):
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
        print self._lModels
        for iLabel in ['tt','st','wlnu']:
            self._lModels[iLabel] = self.getModel(iLabel)
            self._lModels[iLabel].Print()

        # self._lModels['Mc']  = r.RooAddPdf("model_mc","model_mc",
        #                                    r.RooArgList(self._lHPdfs['tt'],self._lHPdfs['st'],self._lHPdfs['wlnu']),
        #                                    r.RooArgList(self._lNPdfs['st'],self._lNPdfs['wlnu']))
        # so apparatenly there should be 2 type of "ttbar" models: model_bkg_data,model_ttbar_data, and also: model_bkg_TotalMC, model_ttbar_TotalMC and I do not understnd the differences...
        self._lModels['Mc'] = r.RooAddPdf("model_mc","model_mc",r.RooArgList(self._lModels['tt'],self._lModels['st'],self._lModels['wlnu']))
        #getattr(self._lW,'import')(self._lModels['Mc'],r.RooFit.RecycleConflictNodes())
        
        self._lModels['Data'] = r.RooAddPdf("model_data","model_data",r.RooArgList(self._lModels['tt'],self._lModels['st'],self._lModels['wlnu']))
        print 'model data'
        self._lModels['Data'].Print()
        

        ## missing constraints? which ones?
        # lDataConstraints = r.RooArgSet("lConstraints_data")
        # for iC in range(self._lConstraints['data']):
        #     lDataConstraints.add(self._lW.pdf(lConstraints[iC]))
        #     lDataConstraints.Print()

        # fit to data
        pRooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE)) #r.ExternalConstraints(lDataConstraints))
        RooFitResult_Data = self._lModels['Data'].fitTo(lSData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
        pRooFitResult_Data.Print()
        draw(lVar,lSData,[self._lModels['Data']],'data_only') # make this another frame? with chi2? and include pull?

        # fit to mc
        pRooFitResult_Mc   = self._lModels['Mc'].fitTo(lSMc,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
        pRooFitResult_Mc   = self._lModels['Mc'].fitTo(lSMc,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE))
        pRooFitResult_Mc.Print()
        draw(lVar,lSMc,[self._lModels['Mc']],'mc_only')

        # draw data and mc
        drawDataMc(lVar,lSData,[self._lModels['Data']],lSMc,[self._lModels['Mc']],"data_mc")

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

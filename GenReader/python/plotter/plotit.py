import ROOT as r 
from inputs import inputs
import psutil
p = psutil.Process()
import copy
lumi=136

def getNormalizedHistoAtPoint(sample, hist, point):
    tf=r.TFile.Open(sample['file'])
    histo=tf.Get(hist)
    if not histo: 
        raise RuntimeError("Histo %s does not exist in %s"%(hist, sample['file']))
    sm_norm=tf.Get("EFTGenHistsWithCuts/h_eventsumEFT").GetBinContent(1, r.WCPoint())

    outhist= histo.Clone()
    outhist.SetDirectory(0)
    for bin in range(1, outhist.GetXaxis().GetNbins()+1):
        fit=histo.GetBinFit( bin )
        
        outhist.SetBinContent( bin, sample['xsec']*lumi * fit.evalPoint( point)  / sm_norm ) 
        outhist.SetBinError  ( bin, sample['xsec']*lumi * fit.evalPointError( point)  / sm_norm ) 
    outhist.SetFillColor(sample['color'])

    tf.Close()
    return outhist


def getStackForPoint( point , histName):
    stack = r.THStack()
    histos=[]
    for sample in inputs:
        hist=getNormalizedHistoAtPoint(sample, histName, point)
        stack.Add(hist)
        histos.append(hist)
    return stack, histos

def getCanvasForRatio():
    plotformat = (1200,600)
    height = plotformat[1]+150 
    c1 = r.TCanvas("canvas", '', plotformat[0], height)
    c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), (plotformat[1]+150 + (plotformat[1]+150 - c1.GetWh())));
    p1 = r.TPad("pad1","pad1",0,0.30,1,1);
    p1.SetTopMargin(p1.GetTopMargin()*1.1);
    p1.SetBottomMargin(0);
    p1.Draw();
    p2 = r.TPad("pad2","pad2",0,0,1,0.30);
    p2.SetTopMargin(0 );
    p2.SetBottomMargin(0.3);
    p2.SetFillStyle(0);
    p2.Draw();
    return (c1, p1,p2)

def plotHistoForPoint( operator, distro):

    smpoint = r.WCPoint()
    stack, hists=getStackForPoint(smpoint, 'EFTGenHistsWithCuts/' + distro )

    bsmpoint=r.WCPoint()
    if operator != "":
        bsmpoint.setStrength(operator,1)
    stack_bsm, hists_bsm=getStackForPoint(bsmpoint, 'EFTGenHistsWithCuts/' + distro)
    if distro =='top19001_cat':
        for h in hists+hists_bsm:
            for idx_label, label in enumerate(["2lss (+)", "2lss (-)", "3l1b (+)", "3l1b (-)", "3l2b (+)", "3l2b (-)", "SFZ1b", "SFZ2b", "4l"]):
                h.GetXaxis().SetBinLabel(idx_label+1, label)
    h.GetXaxis().SetLabelSize(0.2)

    (c1,p1,p2)=getCanvasForRatio()
    c1.cd()
    stack_bsm.Draw("hist")
    c1.SaveAs("%s_%s_bsm.pdf"%(operator,distro))
    c1.SaveAs("%s_%s_bsm.png"%(operator,distro))
    sm=stack.GetStack().Last()
    ratio=copy.deepcopy(stack_bsm.GetStack().Last())
    ratio.Divide(sm)
    ratio.GetYaxis().SetLabelSize(0.1)
    #for i in range(ratio.GetNbinsX()):
    #    print ratio.GetBinContent(i+1)
        
    (c1,p1,p2) = getCanvasForRatio()
    
    p1.cd()
    stack_bsm.Draw("hist")
    stack_bsm.GetHistogram().GetYaxis().SetTitle("Events")
    p2.cd()
    ratio.GetYaxis().SetRangeUser(0,2)
    ratio.Draw()
    c1.SaveAs('%s_%s.pdf'%(distro,operator))
    c1.SaveAs('%s_%s.png'%(distro,operator))

    for h in hists    : h.IsA().Destructor( h ) 
    for h in hists_bsm: h.IsA().Destructor( h ) 

    
if __name__ == "__main__":
    r.gROOT.SetBatch(True)
    r.gROOT.ProcessLine(".x tdrstyle.cc")
    r.gStyle.SetOptStat(0)
    r.gStyle.SetOptTitle(0)

    for k in ["cQq13","cQq83","cQq11","cQq81","ctq1","ctq8"]:
        for distro in ['top19001_cat']:
            plotHistoForPoint(k, distro)
            

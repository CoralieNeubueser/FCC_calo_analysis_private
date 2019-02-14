# Setup ROOT
import ROOT
from ROOT import TH1F, TLegend, TPad, gPad, TCanvas, SetOwnership, TPaveText, TLine

def draw_1histogram( histo, x_axisName, y_axisName ):
   histo.GetXaxis().SetTitle(x_axisName)
   if (y_axisName==""):
      histo.GetYaxis().SetTitle("Entries/per bin")
   else:
      histo.GetYaxis().SetTitle(y_axisName)
   maximum = 1.2*histo.GetMaximum()
   histo.SetMaximum(maximum)
   histo.Draw()
   gPad.Update()
   return

def draw_2histograms( histo1, histo2, x_axisName, y_axisName, leg1Name, leg2Name ):
   histo1.GetXaxis().SetTitle(x_axisName)
   if (y_axisName==""):
      histo1.GetYaxis().SetTitle("Entries/per bin")
   else:
      histo1.GetYaxis().SetTitle(y_axisName)

   histo2.SetLineColor(2)
   
   if (leg1Name!=""):
      histo1.SetStats(0)
      histo2.SetStats(0)
      legend=TLegend(0.45,0.78,0.9,0.9)
      legend.AddEntry(histo1, leg1Name, "l")
      legend.AddEntry(histo2, leg2Name, "l")
      legend.SetTextSize(0.05)
      SetOwnership( legend, 0 )

   maximum1 = 1.2*histo1.GetMaximum()
   maximum2 = 1.2*histo2.GetMaximum()
   maximum = maximum1
   if (maximum2>maximum):
      maximum = maximum2
   histo1.SetMaximum(maximum)
   histo1.Draw()
   histo2.Draw("same")
   if (leg1Name!=""):
      legend.Draw("same")

   gPad.Update()

   return


def draw_1histogram_normalized( histo, x_axisName, y_axisName ):
   hsto.Sumw2()
   histo.Scale(1./histo.GetIngetral())
   draw_1histogram(histo, x_axisName, y_axisName)
   return

def draw_2histograms_normalized( histo1, histo2, x_axisName,y_axisName, leg1Name, leg2Name ):
   histo1.Sumw2()
   histo2.Sumw2()
   histo1.Scale(1./histo1.GetIntegral())
   histo2.Scale(1./histo2.GetIntegral())
   draw_2histograms( histo1, histo2, x_axisName, y_axisName, leg1Name, leg2Name)
   return

def draw_hist2d(hist, titleX = "", titleY = "", titleZ = "", title = ""):
    hist.Draw('colz')
    if titleX:
       hist.GetXaxis().SetTitle(titleX)
    if titleY:
       hist.GetYaxis().SetTitle(titleY)
    if titleZ:
       hist.GetZaxis().SetTitle(titleZ)
    if title:
       hist.SetTitle(title)
       

def draw_text(lines, coordinates = [0.1,0.8,0.5,0.9], colour = 36):
   text = TPaveText(coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    coordinates[3],"brNDC")
   text.SetFillColorAlpha(0,1)
   for line in lines:
      text.AddText("#color["+str(colour)+"]{"+line+"}")
      print(line)
      text.Draw()
      ROOT.SetOwnership(text,False)

def draw_rectangle(start = [0,0], end = [1,1], colour = 2, width = 2):
   lines = []
   lines.append(TLine(start[0],start[1],end[0],start[1]))
   lines.append(TLine(end[0],start[1],end[0],end[1]))
   lines.append(TLine(start[0],end[1],end[0],end[1]))
   lines.append(TLine(start[0],start[1],start[0],end[1]))
   for line in lines:
      line.SetLineColor(colour)
      line.SetLineWidth(width)
      line.Draw()
      ROOT.SetOwnership(line,False)

def draw_resoLin(histReso, histLin):
   pad1 = TPad("pad1","pad1",0,0,1,0.66)
   pad2 = TPad("pad2","pad2",0,0.66,1,1)
   pad2.SetBottomMargin(0.01)
   pad1.SetBorderMode(0)
   pad1.SetTopMargin(0.01)
   pad1.SetBottomMargin(0.15)
   pad2.SetBorderMode(0)
   pad1.SetTickx(1)
   pad2.SetTickx(1)
   pad1.SetTicky(1)
   pad2.SetTicky(1)
   pad1.Draw()
   pad2.Draw()
   
   pad1.cd()
   pad1.SetLogx()
   histReso.Draw("AP")
   histReso.GetXaxis().SetTitle("E [GeV]")
   histReso.GetYaxis().SetTitle("#sigma_{E_{tot}}/#LT E_{tot}#GT")
   histReso.GetXaxis().SetLimits(10,11000)
   histReso.SetMinimum(0.)
   histReso.SetMaximum(0.15)
   histReso.Draw("AP")
 
   pad2.cd()
   pad2.SetLogx()
   pad2.SetGridy()
   histLin.GetYaxis().SetLabelSize(0.09)
   histLin.GetYaxis().SetTitleOffset(0.9)
   histLin.GetYaxis().SetTickLength(0.04)
   histLin.GetYaxis().SetTitleSize(0.09)
   histLin.GetXaxis().SetTickLength(0.04)
   histLin.GetXaxis().SetTitleOffset(1.1)
   histLin.GetXaxis().SetLabelSize(0.1)
   histLin.GetXaxis().SetTitleSize(0.09)
   histLin.GetXaxis().SetLimits(10,11000)
   histLin.GetYaxis().SetRangeUser(-0.155,0.155)
   histLin.Draw("AP")

import csv
from math import log

def th2f_to_csv(hist, csv_file, xmin = 1, xmax = -1, ymin = 1, ymax = -1, scale = 1, saveLog=False):
    """Print TH2F bin data to CSV file."""
    if xmax ==-1:
        xmax = hist.GetNbinsX()
    if ymax ==-1:
        ymax = hist.GetNbinsY()
    xaxis, yaxis = hist.GetXaxis(), hist.GetYaxis()
    with open(csv_file, 'w') as f:
        c = csv.writer(f, delimiter=' ', lineterminator='\n')
        yid = int(-0.5 * (ymax - ymin))
        for ybin in xrange(ymin, ymax+1):
            y_lowedge = yaxis.GetBinLowEdge(ybin)
            xid = int(-0.5 * (ymax - ymin))
            for xbin in xrange(xmin, xmax+1):
                x_lowedge = xaxis.GetBinLowEdge(xbin)
                weight = hist.GetBinContent(xbin, ybin)
                # if weight > 0:
                if saveLog:
                    if weight > 0:
                        c.writerow((x_lowedge, y_lowedge,log( weight / scale) ))
                    else:
                        c.writerow((x_lowedge, y_lowedge, - 20 ))
                else:
                    if weight > 0:
                        c.writerow((xid, yid, weight / scale))
                xid += 1
            yid += 1
        print "SAVED ", (ymax-ymin+1), ' * ', (xmax-xmin+1) , ' = ', (ymax-ymin+1) * (xmax-xmin+1)

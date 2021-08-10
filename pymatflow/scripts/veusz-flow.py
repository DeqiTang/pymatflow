#!/usr/bin/env python

import os
import sys
import copy
import argparse

import veusz.embed as veusz_embed
import time

"""
Note:
    if running in WLS, you need to let Qt not launch the window by:
        export QT_QPA_PLATFORM='offscreen'
"""


def bandDataToVeusz(datafile, gplotfile, vszfile):
    # First you should import the band data to veusz with name 'data' as a nD dataset.

    #  open a new window and return a new Embedded object
    embed = veusz_embed.Embedded('window title', hidden=True)
    
    # read band data file
    # Important need to build another
    embed.ImportFileND(filename=datafile,
            mode="txt", dataset="data")
    
    # extract efermi
    with open(gplotfile, 'r') as fin:
        tmp = fin.readlines()
    for line in tmp:
        if "plot" in line and "using" in line:
            break
    if 1 == line.split("using")[1].split(")")[0].split("$")[1].count("-"):
        efermi = float(line.split("using")[1].split(")")[0].split("$")[1].split("-")[-1])
    elif 2 == line.split("using")[1].split(")")[0].split("$")[1].count("-"):
        efermi = 0 - float(line.split("using")[1].split(")")[0].split("$")[1].count("-")[-1])


    #  make a new page, but adding a page widget to the root widget
    page = embed.Root.Add('page')

    # note: autoAdd=False stops graph automatically adding own axes (used in saved files)
    graph = page.Add('graph', autoadd=False)

    # do the plot of band lines
    data = embed.GetData('data')
    for i in range(data.shape[0]):
        #xy = graph.Add('xy', xData="data[%d, :, 0]" % i, yData="data[%d, :, 1] - %f" % (i, efermi), marker="none")
        xy = graph.Add('xy')
        xy.xData.val ="data[%d, :, 0]" % i
        xy.yData.val ="data[%d, :, 1] - %f" % (i, efermi)
        xy.marker.val ="none"
        xy.PlotLine.width.val = "1.5pt"
    # formatting
    # use widget.childnames to see what properties it has

    x = graph.Add('axis', name='x')
    y = graph.Add('axis', name='y')
    y.direction.val = 'vertical'

    graph.y.label.val = 'Energy (E - E_{f})'
    graph.y.Label.font.val = "Source Han Sans"
    graph.y.TickLabels.font.val = "Source Han Sans"
    graph.y.Label.size.val = "16pt"
    graph.y.TickLabels.size.val = "16pt"
    graph.y.Line.width.val = "3pt"
    graph.x.Line.width.val = "3pt"
    graph.x.TickLabels.hide.val = True
    graph.x.MajorTicks.hide.val = True
    graph.x.MinorTicks.hide.val = True
    graph.y.MajorTicks.width.val = "2pt"
    graph.y.MinorTicks.width.val = "1.5pt"
    graph.y.max.val = 4
    graph.y.min.val  = -4

    # parse gnuplot file for high symmetry kpath
    # find and parse line like:
    #   set xtics('{/symbol G}' 0.000000, 'X' 0.8172, 'M' 1.63562, '{/symbol G}' 2.79211, 'Z' 2.86994, 'R' 3.68714, 'A' 4.50556, 'Z | X' 5.66205, 'R | M' 5.73987, 'A' 5.817703)

    with open(gplotfile, 'r') as fin:
        gp = fin.readlines()
    for i, line in enumerate(gp):
        if "set xtics" in line: # assuming 'set xtics' is in only one line
            break        
    pairs = line.split("\n")[0].replace("set xtics", "").replace("{/symbol G}", "\Gamma").replace("(", "", 1).replace(")", "",-1).split(",")
    klabel = []
    for item in pairs:
        klabel.append({
            "x": float(item.split()[-1]),
            "label": "".join(item.replace("'", "").split()[:-1])
        })

    ## add high symmetry vertical line and label
    for item in klabel:
        line = graph.Add("line")
        line.mode.val = "length-angle"
        line.positioning.val = "axes"
        line.angle.val = [90]
        line.xPos.val = [item["x"]]
        line.yPos.val = [graph.y.max.val]
        line.length.val = [1]
        line.Line.style.val = "dashed"
        line.Line.width.val = "1pt"
        label = graph.Add("label")
        label.Text.font.val = "Source Han Sans"
        label.Text.size.val = "16pt"
        label.positioning.val = "axes"
        label.label.val = item["label"]
        label.xPos.val = [item["x"]]
        label.yPos.val = [graph.y.min.val-0.5]

    graph.x.autoExtend.val = False
    graph.x.autoExtendZero.val = False

    # save to veusz project file
    embed.Save(vszfile)
    embed.Export("%s.png" % vszfile, color=True, page=-1, dpi=300, antialias=True, quality=85, backcolor='#ffffff00', pdfdpi=150, svgtextastext=False)



def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, nargs=2, required=True,
        help="input data file and gnuplot file, -i TOTAL BANDDATA GNUPLOTFILE. if running in WLS, you need to let Qt not launch the window by: export QT_QPA_PLATFORM=\'offscreen\'")

    parser.add_argument("--output-vsz", type=str, default="./auto-output.vsz",
        help="save to an veusz project file")

    args = parser.parse_args()
    #
    #

    bandDataToVeusz(
        datafile=args.input[0],
        gplotfile=args.input[1],
        vszfile=args.output_vsz)

if __name__ == "__main__":
    main()
#!/usr/bin/env python

import os
import sys
import copy
import argparse

import veusz.embed as veusz_embed 
import time

"""
using veusz in python via 'embed' mode
Compatible with Veusz 3.3.1
Note:
    if running in WLS, you need to let Qt not launch the window by:
        export QT_QPA_PLATFORM='offscreen'

Dependencies:
    pip install --user veusz
Ref:
    * https://veusz.github.io/docs/manual/api.html
    * https://github.com/veusz/veusz/wiki
    * https://github.com/jeremysanders/veusz/wiki/EmbeddingPython
    * https://github.com/veusz/veusz/wiki/ToolsPlugins
"""

os.environ["QT_QPA_PLATFORM"] = "offscreen"

def setAllNodesSameFont(node, font='Source Han Sans'):
    """
    Note:
        Nodes: including WidgetNode, SettingGroupNode, SettingNode
        node can start from Root
    """
    # setting the font of every children_widgets children_settinggroups 
    # of Root to the specified font recursively
    for child in list(node.children):
        if 'font' in child.childnames:
            child.font.val = font
        setAllNodesSameFont(child, font=font)

def bandDataToVeuszEmbed(datafile, gplotfile, vszfile, font="Source Han Sans"):
    """
    Note:
        datafile and gplotfile are files usually generated in post-processing dir by pymatflow
    """
    #os.environ["QT_QPA_PLATFORM"] = "offscreen"
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
        efermi = 0 - float(line.split("using")[1].split(")")[0].split("$")[1].split("-")[-1])


    #  make a new page, but adding a page widget to the root widget
    page = embed.Root.Add('page')
    # if you run python in Veusz GUI application at console windows,
    # there is a default global 'Root' already which handle the current GUI.
    # and the command like ImportFileND are all available.
    # However, we adopt 'embed' mode here, all the above are available under
    # the embed.Embedded() class.

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
        #xy.PlotLine.width.val = "1.5pt"
        xy.PlotLine.width.val = "2pt"
        xy.PlotLine.color.val = "black"
    # formatting
    # use widget.childnames to see what properties it has

    x = graph.Add('axis', name='x')
    y = graph.Add('axis', name='y')
    y.direction.val = 'vertical'

    graph.y.label.val = 'Energy (E - E_{f})'
    graph.y.Label.font.val = font
    graph.y.TickLabels.font.val = font
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
    graph.leftMargin.val = "2cm"
    graph.rightMargin.val = "0.5cm"
    graph.aspect.val = 1 # default is 'auto'

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
        line.positioning.val = "relative"
        line.angle.val = [90]
        line.xPos.val = [item["x"] / klabel[-1]["x"]]
        line.yPos.val = [1]
        line.length.val = [1]
        line.Line.style.val = "dashed"
        line.Line.width.val = "1pt"
        label = graph.Add("label")
        label.Text.font.val = font
        label.Text.size.val = "16pt"
        label.positioning.val = "relative"
        label.label.val = item["label"]
        # minus by shiftLabelX to set the label aligned with the verticalline by vertical center
        # because by dfault the label will set its left bound to label.xPos.val
        shiftLabelX = float(label.Text.size.val.split("pt")[0]) / 15 * 0.015
        label.xPos.val = [item["x"] / klabel[-1]["x"] - shiftLabelX] 
        label.yPos.val = [0 - 0.05] # a little bit lower than the bottom x axis
    # and Fermi horizontal line
    line = graph.Add("line")
    line.mode.val = "length-angle"
    line.positioning.val = "relative"
    line.angle.val = [0]
    line.xPos.val = [0]
    line.yPos.val = [abs(0-graph.y.min.val) / (graph.y.max.val - graph.y.min.val)]
    line.length.val = [1]
    line.Line.style.val = "dashed"
    line.Line.width.val = "0.5pt"

    graph.x.autoExtend.val = False
    graph.x.autoExtendZero.val = False

    # set all font to font
    setAllNodesSameFont(embed.Root, font=font)

    # save to veusz project file
    embed.Save(vszfile)
    #embed.Export("%s.png" % os.path.join(os.path.dirname(vszfile), os.path.basename(vszfile)), color=True, page=-1, dpi=300, antialias=True, quality=100, backcolor='#ffffff00', pdfdpi=150, svgtextastext=False)
    embed.Export("%s.png" % vszfile, color=True, page=-1, dpi=300, antialias=True, quality=100, backcolor='#ffffff00', pdfdpi=150, svgtextastext=False)

    #del os.environ["QT_QPA_PLATFORM"]

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, nargs=2, required=True,
        help="input data file and gnuplot file usually found in post-processing dir from pymatflow, -i TOTAL BANDDATA GNUPLOTFILE. if running in WLS, you need to let Qt not launch the window by: export QT_QPA_PLATFORM=\'offscreen\'")

    parser.add_argument("--output-vsz", type=str, default="./auto-output.vsz",
        help="save to an veusz project file")

    parser.add_argument("--font", type=str, default="Arial",
        #choices=["Source Han Sans", "Aria", "Times New Roman", "Noto Sans CJK SC"],
        help="set the font")
    args = parser.parse_args()
    #
    #

    bandDataToVeuszEmbed(
        datafile=args.input[0],
        gplotfile=args.input[1],
        vszfile=args.output_vsz,
        font=args.font
    )

if __name__ == "__main__":
    main()
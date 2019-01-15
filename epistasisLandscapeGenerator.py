#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 15:54:16 2019

@author: Caleb Renshaw

epistasisLandscapeGenerator.py - Generates a diagram displaying epistatic relations between alleles.

Takes a .csv file where the first column is a strain identifier, and the next n columns are the alleles.

filePath is the relative or absolute path to the .csv file containing the data.

fitCol is the .csv column number (0-indexed) with the relative fitness data to plot.

pathSigs is a list of strings of length n representing the genotype signature for each pathway to plot,
  use '0' and '1' to represent ancestral and evolved alleles, '*' to represent either, 
  and '2' to represent neither e.g.: '****2' '****1' '2***1'; use 0, 1, & 2 in the .csv file as well.

pathNames is a list of names for each pathway; must be same length as pathSigs.

Ex_Usage: $ python3 epistasisLandscapeGenerator.py <filePath> <fitCol> --pathSigs 'sig_1' 'sig_2' --pathNames 'name_1' 'name_2'

Note: this script should ideally not be run through Anaconda on Ubuntu Linux, as modern fonts accessed with tkinter will not be rendered.
Use the system's python3 instead.
"""

# if any of the following modules are missing, run from terminal:
# pip install [package]  (or  pip3 install [package]  if that doesn't work)

import tkinter as tk  # tkinter should come with python, but if not, use sudo apt-get install python3-tk
import argparse, re, datetime, os # likewise should already come with python, possible exception of re
import pandas as pd
import pyscreenshot as ImageGrab
#from PIL import ImageGrab  <- use this for Windows instead of pyscreenshot

CLI = argparse.ArgumentParser()
CLI.add_argument('dataFilename', nargs=1, type=str)
CLI.add_argument('fitCol', nargs=1, type=int)
CLI.add_argument('--pathSigs', nargs='*', type=str)
CLI.add_argument('--pathNames', nargs='*', type=str)
args = CLI.parse_args()

print('dataFilename: %r' % args.dataFilename)
print('fitCol: %r' % args.fitCol)
print('pathSigs: %r' % args.pathSigs)
print('pathNames: %r' % args.pathNames)

dataFilename = args.dataFilename[0]
fitCol = args.fitCol[0]
pathSigs = args.pathSigs
pathNames = args.pathNames

dataFilename = os.path.abspath(dataFilename)

test = lambda x: len(set(map(len, x))) == 1     # check that all pathSigs elements are the same length
assert(test(pathSigs)), 'Elements in pathSigs need to be the same length.'
assert(len(pathSigs) == len(pathNames)), 'pathSigs and pathNames need to be the same length.'

df = pd.read_csv(dataFilename)

numPaths = len(pathSigs)    # number of pathways to map
numGenes = len(pathSigs[0]) # number of genes
genes = []  # list of genes to map
cols = [0]  # list of dataframe columns of interest
for i in range(1, numGenes+1):
    genes.append(list(df)[i])
    cols.append(i)
cols.append(fitCol)

# default color pallet for first 10 genes and first 6 pathways
colors = ['red', 'magenta2', 'forest green', 'blue', 'black', 'gold', 'lime green', 'saddle brown', 'dark violet', 'cornflower blue']
pathColors = ['orange', 'purple', 'cyan', 'khaki', 'olive drab', 'royal blue']
geneCols = {}   
pathCols = {}

for i in range(len(genes)):
    try:
        geneCols[genes[i]] = colors[i]
    except KeyError:    # snark
        print('Not enough colors assigned by default (10) for genes provided.')
        print('Add more colors to the default pallet "colors" to continue.')
        print('But really, you should reconsider mapping that many genes to begin with...')

for i in range(numPaths):
    try:
        pathCols[pathNames[i]] = pathColors[i]
    except KeyError:    # snark
        print('Not enough colors assigned by default (6) for pathways provided.')
        print('Add more colors to the default pallet "pathColors" to continue.')
        print('But really, you should reconsider mapping that many pathways to begin with...')

# define binomial theorem for determining nodes per pathway
def binomial(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

nodeCount = []  # list of nodes per layer in a pathway of n variable alleles
def nodes(n):
    nodeCount.clear()
    for i in range(n+1):
        nodeCount.append(binomial(n, i))
    return nodeCount

strains = []    # Generate list of single-string strain genotype IDs
for i in range(len(df.index)):
    strains.append(''.join(map(str, list(df.iloc[i,1:(numGenes+1)]))))
assert(len(strains) == len(set(strains))), 'Sourcefile cannot contain duplicate genotypes.'

pstrains = []   # Generate identical list, but with '-' substituted for '2'; used for display
for i in range(len(strains)):
    pstrains.append(re.sub('2', '-', strains[i]))

# Separate main dataframe into subsets per pathway
paths = {'indices': {},     # indicies of the main dataframe per pathway
         'dfs': {},         # dataframe subsets per pathway
         'coords': {},      # coordinates for plotting the pathway nodes
         'cxns': {}}        # connections between nodes, along with +/- indication
varAlleles = []       # number of variable alleles per pathway
mutCount = []         # number of evolved/mutant alleles per strain
pathSpanNodes = []    # list of max nodes to span horizontally per pathway

for i in range(numPaths):
    # turn pseudo regexes from pathSigs to valid regex
    pathSig = re.sub('\*', '(0|1)', pathSigs[i])
    pathRegex = re.compile(pathSig)
    # build data structure
    paths['indices'][i] = [j for j, x in enumerate(strains) if re.search(pathRegex, x)]
    paths['dfs'][i] = df.iloc[paths['indices'][i], cols]
    paths['dfs'][i].rename(columns={paths['dfs'][i].columns[numGenes+1]: 'Fitness'}, inplace=True)
    paths['dfs'][i].iloc[:,numGenes+1] = pd.to_numeric(paths['dfs'][i].iloc[:,numGenes+1])
    paths['dfs'][i]['Genotype'] = [strains[j] for j in paths['indices'][i]]
    paths['dfs'][i]['pGenotype'] = [pstrains[j] for j in paths['indices'][i]]
    paths['coords'][i] = {}
    paths['cxns'][i] = {}
  
    var = 0
    const = []
    for j in range(1, numGenes+1):  # identify constant and variable alleles per path
        if len(set(paths['dfs'][i].iloc[:,j])) != 1:
            var += 1
        else:
            const.append(j)
    
    mutCount.clear()
    for j in range(len(paths['dfs'][i])):
        try:    # count mutant alleles per strain
            mutCount.append(paths['dfs'][i].iloc[j,1:(numGenes+1)].value_counts()[1])
        except KeyError:    # return 0 if no mutants are found
            mutCount.append(0)

    paths['dfs'][i]['mutCount'] = mutCount
    paths['dfs'][i].sort_values(by=['mutCount', 'Genotype'], ascending=[True, False], inplace=True)
    paths['dfs'][i].index = range(len(paths['dfs'][i]))
        # sort the dataframes first by mutation count (low-high), then by genotype (high-low)
        # for printing pathways bottom-top then left-right
    varAlleles.append(var)
    pathSpanNodes.append(max(nodes(varAlleles[i])))

nodesVert = max(varAlleles) + 1     # total number of nodes to span vertically

pathBases = [list(paths['dfs'][x].loc[0,genes]) for x in range(numPaths)]
pathLinks = {}      # connections between pathway bases

for i in range(numPaths - 1):   # assumes pathways will be connected only to immediate neighbors
    for j in range(numGenes):
        if (pathBases[i][j] == 2) & (pathBases[i+1][j] in [0, 1]):      # addition of a gene
            pathLinks[i] = ['+', list(paths['dfs'][i])[j+1]]
        elif (pathBases[i][j] in [0, 1]) & (pathBases[i+1][j] == 2):    # deletion of a gene
            pathLinks[i] = ['-', list(paths['dfs'][i])[j+1]]

allDifs = []    # list of fitness differences between adjacent nodes, to normalize
for i in range(numPaths):   # establish connections between nodes per path
    paths['cxns'][i] = {}
    
    for j in range(len(paths['dfs'][i])):   # for each ancestral strain in a path
        paths['cxns'][i][j] = {}
        hits = []       # all strains precisely one mutation more evolved than the ancestor
        alleles = []    # which alleles the ancestor may evolve
    
        for k in range(1, numGenes+1):
            if paths['dfs'][i].iloc[j, k] == 0:     # if an allele is ancestral at a locus
                allele = list(paths['dfs'][i])[k]   # identify the allele
                paths['dfs'][i].iloc[j, k] = 1      # briefly change to evolved genotype
                seq = ''.join(map(str, list(paths['dfs'][i].iloc[j,1:(numGenes+1)])))
                hit = paths['dfs'][i].index[paths['dfs'][i]['Genotype'] == seq].tolist()    # find strain evolved at locus
                assert(len(hit) == 1), 'Something went wrong, please check genotypes for errors' # should only be one match
                paths['dfs'][i].iloc[j, k] = 0      # return allele to normal
                hits += hit
                alleles.append(allele)

        paths['cxns'][i][j]['hits'] = hits
        paths['cxns'][i][j]['alleles'] = alleles
        paths['cxns'][i][j]['fits'] = [round(x , 3) for x in list(paths['dfs'][i].loc[hits, 'Fitness'])]    # list of fitness of strains evolved from ancestor
        paths['cxns'][i][j]['difs'] = [round(x - paths['dfs'][i].loc[j, 'Fitness'], 3) for x in paths['cxns'][i][j]['fits']]    # list of fitness differences
        allDifs += paths['cxns'][i][j]['difs']

normDifs = {}   # normalized fitness differences for plotting line thicknesses
for i in allDifs:
    normDifs[i] = (i - min(allDifs)) / (max(allDifs) - min(allDifs))
normDifs['zero'] = (-min(allDifs)) / (max(allDifs) - min(allDifs))  # normalized 0 value, above is benefit, below is detriment

#####
# data is prepared: now to build objects and determine coordinates
#####

tbw = 60    # text box width
tbh = 40    # text box height
canw = tbw * sum(pathSpanNodes) * 2  # canvas width
canh = tbh * nodesVert * 3           # canvas height

span_x = []     # top-left x-coordinates for max number of nodes to evenly span horizontally
for i in range(sum(pathSpanNodes)):
    span_x.append(int((tbw * i * 2) + (tbw / 2)))   # horizontal distance between nodes is as wide as a node

pathSpan_x = {}     # top-left x-coordinates of longest node row per pathway
start = end = 0
for i in range(len(pathSpanNodes)):
    end += pathSpanNodes[i]
    pathSpan_x[i] = span_x[start:end]
    start += pathSpanNodes[i]
    
pathMid_x = {}  # top-left x-coordinate of middle node
for i in range(len(pathSpan_x)):
    pathMid_x[i] = int((pathSpan_x[i][-1] + pathSpan_x[i][0]) / 2)

span_y = {}     # heights of nodes to span vertically
path_x = {}     # left x-coordinates of each node per pathway, arranged by vertical level
for i in range(numPaths):
    l = 0   # unique index for each node in a path
    span_y[i] = []
    path_x[i] = {}

    for j in range(len(nodes(varAlleles[i]))):  # vertical levels
        span_y[i].append(canh - (tbh + (tbh * j * 2.5))) # vertical distance between nodes is a node height-and-a-half
        path_x[i][j] = []

        if nodes(varAlleles[i])[j] % 2 == 1:    # if j is a layer with odd nodes
            flanks = int((nodes(varAlleles[i])[j] - 1) / 2) 
            if flanks == 0:
                path_x[i][j] = [pathMid_x[i]]
            else:
                for k in range(1, flanks + 1): # node x-coords are node-width multiples from middle node
                    lFlank = [pathMid_x[i] - (tbw * k * 2)]
                    rFlank = [pathMid_x[i] + (tbw * k * 2)]
                path_x[i][j] = lFlank + [pathMid_x[i]] + rFlank

        else:   # if j is a layer with even nodes
            flanks = int((pathSpanNodes[i] - nodes(varAlleles[i])[j]) / 2)
            if flanks == 0:     # layer is widest layer
                path_x[i][j] = pathSpan_x[i]
            else:   # x-coords are equal to the middle values of widest layer
               path_x[i][j] = pathSpan_x[i][flanks:-flanks]

        for k in range(nodes(varAlleles[i])[j]):    # nodes within levels
            paths['coords'][i][l] = {'x1': path_x[i][j][k], 'y1': span_y[i][j], 'x2': path_x[i][j][k] + tbw, 'y2': span_y[i][j] - tbh,
                 'xc': ((path_x[i][j][k] * 2) + tbw) / 2, 'yc': ((span_y[i][j] * 2) - tbh) / 2}  
                 # coordinates of each node corner and central x- and y-coordinate
            l += 1

min_y = canh    # y-coord of uppermost node; (0,0) is top-left corner of Canvas
for i in span_y:
    if min(span_y[i]) < min_y:
        min_y = min(span_y[i])
    
    
cornKeys = ['x1', 'y1', 'x2', 'y2']     # rectangle corner (top-left, bottom-right) coordinates
cenTKeys = ['xc', 'y2']                 # rectangle center-top coordinate
cenBKeys = ['xc', 'y1']                 # rectangle center-bottom coordinate
cenLKeys = ['x1', 'yc']                 # rectangle center-left coordinate
cenRKeys = ['x2', 'yc']                 # rectangle center-right coordinate


# Create a canvas to display the landscape maps
class App(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        self.cv = tk.Canvas(self, width=canw, height=canh)
        
        for i in range(numPaths):
            for j in range(len(paths['coords'][i])):    # plot a rectangle containing a genotype and fitness for each node
                for k in range(len(paths['cxns'][i][j]['hits'])):   # draw lines between adjacent nodes
                    lineStart = tuple([paths['coords'][i][j][x] for x in cenTKeys])
                    e = paths['cxns'][i][j]['hits'][k]
                    lineEnd = tuple([paths['coords'][i][e][x] for x in cenBKeys])
                    gColor = geneCols[paths['cxns'][i][j]['alleles'][k]] # line color based on allele
                    if normDifs[paths['cxns'][i][j]['difs'][k]] <= normDifs['zero']:
                        dashPat = (4,4)     # line is dashed if mutation is detrimental
                        lineWidth = -10 * (normDifs[paths['cxns'][i][j]['difs'][k]] - normDifs['zero'])
                    else:
                        dashPat = None      # line is solid if mutation is beneficial
                        lineWidth = 10 * (normDifs[paths['cxns'][i][j]['difs'][k]] - normDifs['zero'])
                            # line thickness is a function of fitness change
                    self.cv.create_line(lineStart, lineEnd, dash=dashPat, fill=gColor, width=lineWidth, capstyle=tk.ROUND)

                recCoords = tuple([paths['coords'][i][j][x] for x in cornKeys])
                text_x = (recCoords[0] + recCoords[2]) / 2
                text_y = (recCoords[1] + recCoords[3]) / 2 
                geno = paths['dfs'][i].loc[j, 'pGenotype']  # Genotype to print in node
                fit = round(paths['dfs'][i].loc[j, 'Fitness'], 3)   # Fitness to print in node
                self.cv.create_rectangle(recCoords, fill='white')
                self.cv.create_text(text_x, text_y - 10, text=geno, font='arial 15 bold')
                self.cv.create_text(text_x, text_y + 10, text=fit, font='arial 13')
                # apparently anaconda has a hard time dealing with modern fonts accessed by tkinter
                # if this script is run through anaconda, the text will likely be ugly

            pathLabel = pathNames[i]    # place a pathway name over each pathway
            labelColor = pathCols[pathNames[i]]
            text_x = pathMid_x[i] + (tbw / 2)
            text_y = min_y - (tbh + 30)
        
            self.cv.create_text(text_x, text_y, text=pathLabel, font='arial 24 italic', fill=labelColor)
            
        for i in list(pathLinks.keys()):    # lines connecting pathways
            lineStart = tuple([paths['coords'][i][0][x] for x in cenRKeys])
            lineEnd = tuple([paths['coords'][i+1][0][x] for x in cenLKeys])
            gColor = geneCols[pathLinks[i][1]]
            alleleChange = ''.join(pathLinks[i]) # label the gene being added or removed
            text_x = (lineStart[0] + lineEnd[0]) / 2
            text_y = lineStart[1] - 20
    
            self.cv.create_line(lineStart, lineEnd, fill=gColor, width=3)
            self.cv.create_text(text_x, text_y, text=alleleChange, font='arial 20 italic', fill=gColor)
        
        self.cv.create_text(20, canh - (tbh / 2), text='Key:', font='arial')
        for i in range(numGenes):
            gName = genes[i]
            gColor = geneCols[genes[i]]
            text_x = (tbw * (i+1))
            text_y = canh - (tbh / 2)
            self.cv.create_text(text_x, text_y, text=gName, fill=gColor, font='arial')
        for i in range(numGenes-1):
            slash_x = (tbw * (i+1)) + (tbw / 2)
            slash_y = canh - (tbh / 2)
            self.cv.create_text(slash_x, slash_y, text='/', font='arial')
        
        self.cv.grid(row=0, column=0, columnspan=3, sticky='nsew')
        
        # add a button to save a screenshot of the Canvas to the current working directory
        self.snapsave = tk.Button(self, text='Save image to current working directory.', command=self._snapsaveCanvas)
        self.snapsave.grid(row=2, column=0, columnspan=1, sticky='nsew')
        
    def _snapsaveCanvas(self):
        canvas = self._canvas()
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        stampname = ('landscape-' + timestamp + '.png')  # add a timestamp marking when a screenshot was saved
        self.grabcanvas = ImageGrab.grab(bbox=canvas).save(stampname)
        print('Screenshot of canvas image saved as ' + stampname + '.')
        
    def _canvas(self):
        x = self.cv.winfo_rootx() + self.cv.winfo_x()
        y = self.cv.winfo_rooty() + self.cv.winfo_y()
        x1 = x + self.cv.winfo_width()
        y1 = y + self.cv.winfo_height()
        box = (x, y, x1, y1,)
        return box
    
if __name__ == '__main__':
    root = tk.Tk()
    root.title('Fitness Landscapes')
    app = App(root)
    app.grid(row=0, column=0, sticky='nsew')
    
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    
    app.rowconfigure(0, weight=10)
    app.rowconfigure(1, weight=1)
    app.columnconfigure(0, weight=1)
    app.columnconfigure(1, weight=1)
    app.columnconfigure(2, weight=1)
    
    app.mainloop()

# while the Canvas window is open, the terminal which launched the script will be occupied

import ipyrad.analysis as ipa
import toyplot

# the path to your .snps.hdf5 database file
data = "PGL_PHO.snps.hdf5"

# group individuals into populations
imap = {
	"AUS": ["14716_TAS", "14717_TAS", "14719_TAS", "14720_TAS", "14721_TAS", "14751_VIC", "14753_VIC", "14763_VIC", "14773_VIC", "14794_VIC", "14795_VIC", "14802_VIC", "14817_VIC"],
   "NZ": ["15941_NZN", "15951_NZN", "15953_NZN", "15966_NZN", "15976_NZN", "15992_NZN", "16051_NZS", "16066_NZS", "16090_NZS", "16109_NZS", "16136_NZS", "17207_NZS"],
	"CHI": ["15532_CHI", "15551_CHI", "15553_CHI", "15561_CHI", "15564_CHI", "15566_CHI", "15569_CHI", "15572_CHI", "15573_CHI", "15574_CHI", "15596_CHI", "15599_CHI", "15602_CHI"],
	"PHO": ["15910_NZN", "15911_NZN", "15916_NZN", "16000_NZS", "16002_NZS", "16003_NZS", "16009_NZS", "16010_NZS", "16012_NZS", "16085_NZS", "16138_PHO", "16139_PHO", "16140_PHO", "16141_PHO"],
	
}

# require that 50% of samples have data in each group
minmap = {i: 0.5 for i in imap}

# init analysis object with input data and (optional) parameter options
struct = ipa.structure(
    name="test_large",
    data=data,
    imap=imap,
    minmap=minmap,
    mincov=0.9,
)

########################################################
Samples: 52
Sites before filtering: 50984
Filtered (indels): 0
Filtered (bi-allel): 1200
Filtered (mincov): 42983
Filtered (minmap): 36440
Filtered (subsample invariant): 361
Filtered (minor allele frequency): 0
Filtered (combined): 43334
Sites after filtering: 7879
Sites containing missing values: 6271 (79.59%)
Missing values in SNP matrix: 17537 (4.28%)
SNPs (total): 7879
SNPs (unlinked): 2653
########################################################

## Run STRUCTURE and plot results

struct.mainparams.burnin = 5000
struct.mainparams.numreps = 10000

struct.run(nreps=3, kpop=[2, 3, 4, 5], auto=True)

## Analyze results: Choosing K

etable = struct.get_evanno_table([2, 3, 4, 5])
etable

####################################################################################
   Nreps        lnPK       lnPPK   deltaK  estLnProbMean  estLnProbStdev
2      3       0.000       0.000    0.000     -50299.933         133.278
3      3    6833.133    1121.967    1.332     -43466.800         842.199
4      3    5711.167  288945.700  119.305     -37755.633        2421.918
5      3 -283234.533       0.000    0.000    -320990.167      499335.817
####################################################################################

# get canvas object and set size
canvas = toyplot.Canvas(width=400, height=300)

# plot the mean log probability of the models in red
axes = canvas.cartesian(ylabel="estLnProbMean")
axes.plot(etable.estLnProbMean * -1, color="darkred", marker="o")
axes.y.spine.style = {"stroke": "darkred"}

# plot delta K with its own scale bar of left side and in blue
axes = axes.share("x", ylabel="deltaK", ymax=etable.deltaK.max() + etable.deltaK.max() * .25)
axes.plot(etable.deltaK, color="steelblue", marker="o");
axes.y.spine.style = {"stroke": "steelblue"}

# set x labels
axes.x.ticks.locator = toyplot.locator.Explicit(range(len(etable.index)), etable.index)
axes.x.label.text = "K (N ancestral populations)"

import toyplot.pdf
toyplot.pdf.render(canvas, "figure1_large.pdf")


## Analyze results: Barplots

k = 2
table = struct.get_clumpp_table(k)

# sort list by columns
table.sort_values(by=list(range(k)), inplace=True)

# or, sort by a list of names (here taken from imap)
import itertools
onames = list(itertools.chain(*imap.values()))
table = table.loc[onames]

# build barplot
canvas = toyplot.Canvas(width=777, height=333)
axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
axes.bars(table)

# add labels to x-axis
ticklabels = [i for i in table.index.tolist()]
axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
axes.x.ticks.labels.angle = -60
axes.x.ticks.show = True
axes.x.ticks.labels.offset = 10
axes.x.ticks.labels.style = {"font-size": "10px"}

import toyplot.pdf
toyplot.pdf.render(canvas, "figure2_k2_large.pdf")

k = 3
table = struct.get_clumpp_table(k)

# sort list by columns
table.sort_values(by=list(range(k)), inplace=True)

# or, sort by a list of names (here taken from imap)
import itertools
onames = list(itertools.chain(*imap.values()))
table = table.loc[onames]

# build barplot
canvas = toyplot.Canvas(width=777, height=333)
axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
axes.bars(table)

# add labels to x-axis
ticklabels = [i for i in table.index.tolist()]
axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
axes.x.ticks.labels.angle = -60
axes.x.ticks.show = True
axes.x.ticks.labels.offset = 10
axes.x.ticks.labels.style = {"font-size": "10px"}

import toyplot.pdf
toyplot.pdf.render(canvas, "figure2_k3_large.pdf")

k = 4
table = struct.get_clumpp_table(k)

# sort list by columns
table.sort_values(by=list(range(k)), inplace=True)

# or, sort by a list of names (here taken from imap)
import itertools
onames = list(itertools.chain(*imap.values()))
table = table.loc[onames]

# build barplot
canvas = toyplot.Canvas(width=777, height=333)
axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
axes.bars(table)

# add labels to x-axis
ticklabels = [i for i in table.index.tolist()]
axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
axes.x.ticks.labels.angle = -60
axes.x.ticks.show = True
axes.x.ticks.labels.offset = 10
axes.x.ticks.labels.style = {"font-size": "10px"}

import toyplot.pdf
toyplot.pdf.render(canvas, "figure2_k4_large.pdf")

k = 5
table = struct.get_clumpp_table(k)

# sort list by columns
table.sort_values(by=list(range(k)), inplace=True)

# or, sort by a list of names (here taken from imap)
import itertools
onames = list(itertools.chain(*imap.values()))
table = table.loc[onames]

# build barplot
canvas = toyplot.Canvas(width=777, height=333)
axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
axes.bars(table)

# add labels to x-axis
ticklabels = [i for i in table.index.tolist()]
axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
axes.x.ticks.labels.angle = -60
axes.x.ticks.show = True
axes.x.ticks.labels.offset = 10
axes.x.ticks.labels.style = {"font-size": "10px"}

import toyplot.pdf
toyplot.pdf.render(canvas, "figure2_k5_large.pdf")
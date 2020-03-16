import os
import sys
from io import StringIO
import pkg_resources
from seqplot.p9tools import geom_seq_x
from plotnine import ggplot,geom_rect, geom_point, aes, stat_smooth,geom_bar, xlim, ylim, facet_wrap, theme_bw,theme_xkcd, geom_line, geom_tile
from plotnine import theme, scale_y_continuous,scale_x_continuous, theme_bw,theme_classic, theme_dark, theme_light, theme_matplotlib, theme_minimal, theme_seaborn, theme_void
import pandas as pd
import numpy as np



def test_geom_seq_x():

	DATA_PATH = pkg_resources.resource_filename('seqplot', 'data/')
	img=open(os.path.join(DATA_PATH,'seq.png'), 'rb').read()
	os.system('mkdir -p test_results')
	df=pd.DataFrame({'x':np.arange(103),'yy':-np.abs(np.sin(np.arange(103)/10.))})
	g=(ggplot(data=df,mapping=aes(x='x', y='yy')) + geom_point(size=0.1) + geom_bar(stat='identity') + geom_seq_x(seqimg=img,x=df.x,y=df.yy,aspect_ratio=0.1) + theme_light()+theme(aspect_ratio=0.1,dpi=300,plot_margin=0))
	g.save('test_results/plot.png')
	size=os.path.getsize('test_results/plot.png')
	assert size>1000, "output png filesize too small, looks that nothing was produced"


#     assert size>10000, "output png filesize too small, looks that nothing was produced"





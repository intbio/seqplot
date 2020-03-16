"""
This is an advance wrapper in seqplot that helps tp work easily with data that is mapped to sequences from pdb structures.
The residues in pdb chains are numbered with residue ids (resid parameter in MDAnalysis and many other programs).
The chains have also their IDs (segid in MDAnalysis terminology)
The first goal is:
given a pdb id (e.g. 1KX5) and a dataframe of data in the form
segid,resid, value
produce plots of params ontop of the sequences.

"""

from plotnine.utils import resolution
from plotnine.doctools import document
from plotnine.geoms.geom_tile import geom_tile
from plotnine import aes, theme

from pytexshade import ipyshade
from plotnine import ggplot,geom_rect, geom_point, aes, stat_smooth,geom_bar, xlim, ylim, facet_wrap, theme_bw,theme_xkcd, geom_line, geom_tile
from plotnine import facet_wrap, theme, scale_y_continuous,scale_x_continuous, theme_bw,theme_classic, theme_dark, theme_light, theme_matplotlib, theme_minimal, theme_seaborn, theme_void
from .p9tools import geom_seq_x

import pandas as pd
import numpy as np



def plot_prof4pdb(pdbchid='1KX5_A',column='value',data=None,ymin=None):
    """
    Plot profiles on PDB sequences
   a dataframe of data in the form segid,resid, value
    """
    if not isinstance(data,pd.DataFrame):
        data=pd.DataFrame({'segid':['A']*135,'resid':np.arange(1,136),'value':np.sin(np.arange(1,136))})
    
    plot=(ggplot(data=data,mapping=aes(x='resid', y=column))
        + geom_point(size=0.1)+geom_bar(stat='identity')
        + scale_x_continuous(expand=(0,0),name='',breaks=[])
       # + scale_y_continuous(breaks=[0,0.5,1.0])
        + theme_light()+theme(aspect_ratio=0.15,dpi=600,plot_margin=0)) #+ facet_wrap('~ segid',dir='v')

    shaded=ipyshade.shadepdbquick(pdb_chain_id=pdbchid,entrez_email='info@example.com',debug=False,\
                                  force_feature_pos='bottom',ruler='bottom',legend=False,\
                                  feature_types=['SecStr'],show_seq_names=False,show_seq_length=False,density=150)
    
    
    plot = plot + geom_seq_x(seqimg=shaded.img,x=data.resid,y=data[column],ymin=ymin,aspect_ratio=0.15)
    
    
    return plot
        
        
        


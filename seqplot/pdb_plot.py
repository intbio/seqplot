"""
This is an advanced wrapper in seqplot that helps to work easily with data that is mapped to sequences from pdb structures.
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

from pytexshade import ipyshade,shade
from plotnine import ggplot,geom_rect, geom_point, aes, stat_smooth,geom_bar, xlim, ylim, facet_wrap, theme_bw,theme_xkcd, geom_line, geom_tile
from plotnine import facet_wrap, theme, scale_y_continuous,scale_x_continuous, theme_bw,theme_classic, theme_dark, theme_light, theme_matplotlib, theme_minimal, theme_seaborn, theme_void
from .p9tools import geom_seq_x

import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio import SeqIO
from Bio import Entrez
from Bio.Align import MultipleSeqAlignment

import requests
import io

def plot_prof4pdb(pdb_chain_id='1KX5_A',column='value',data=None,ymin=0,draft=False,entrez_email='info@example.com',feature_types=['all'],shading_modes=['charge_functional'],right_overhang_fix=None,debug=False,startnumber=1,cropseq=(None,None)):
    """
    Plot profiles on PDB sequences with annotations
   a dataframe of data in the form segid,resid, value
    """
    if draft:
        dpi=100
    else:
        dpi=300
    if not isinstance(data,pd.DataFrame):
        data=pd.DataFrame({'segid':['A']*135,'resid':np.arange(1,136),'value':np.sin(np.arange(1,136))})
    
    #The problem is that shadepdbquick uses NCBI to get SEQREC, and our data is in resid numbering
    #We need somehow to make correspondence between the two.
    #And ultimately understand where resid=1 is on the SEQREC record, 
    # i.e. what to add to resid, so that resid of seqrec N-terminus will be number 1.
    #To accomplish this we need:
    # 1)real pdb-seq of the chain and its starting resid, also fill gapps with X
    # 2) The seqrec of the chain
    # 3) align pdb-seq to seqrec with no gaps
    # 4) get the length of the N-terminus overhang of seqrec
    # 5) postulate that first resid number in chain is now (overhang_length + 1) (we use 1-based numbering here)
    # The seqrec we can get from NCBI or via Biopython.
    # How to get pdb-seq, looks like Biopython PdbIO is also good,\
    # because we get sequence with Xs at gaps and also the start and end resid numbers.
    # minus is that it does not work with DNA or RNA, but we can make a pull request to biopython with a feature in the future
    if (startnumber!=1):
        print('Warning: stratnumber parameber is known to be buggy, e.g. 0 and -1 values do not work as expected. Check!')
    pdbid=pdb_chain_id.split('_')[0]
    chid=pdb_chain_id.split('_')[1]
    seqrec={}
    pdbseq={}
    h=io.StringIO(requests.get('https://files.rcsb.org/download/%s.pdb'%pdbid).content.decode("utf-8") )
    for record in SeqIO.parse(h,'pdb-seqres'):
        seqrec[record.id.split(':')[1]]=record

    h=io.StringIO(requests.get('https://files.rcsb.org/download/%s.pdb'%pdbid).content.decode("utf-8") )
    for record in SeqIO.parse(h,'pdb-atom'):
        pdbseq[record.id.split(':')[1]]=record
    resid_start=pdbseq[chid].annotations['start']
    alignment = pairwise2.align.globalxs(seqrec[chid].seq[cropseq[0]:cropseq[1]], pdbseq[chid].seq,-10,-1,penalize_end_gaps=False)[0]
#     print(alignment)
        
    #TODO: It's good to add some asserts here to check if the sequence in data corresponds to what we have in seqrec (?)
    
    overhang=len(alignment[1].split(pdbseq[chid].seq[0])[0])
    
    
    datafixed=data
    datafixed['resid']=datafixed['resid']-resid_start+overhang+1
#     print(datafixed)
    #proceed with plotting
    
    #Let's get annotation from genbank indluding secondary structure
    
    
    Entrez.email = entrez_email  # Always tell NCBI who you are
    handle = Entrez.efetch(db="protein", id=pdb_chain_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    msar=MultipleSeqAlignment([record])
    msar[0].id='PDB_'+pdb_chain_id

    msar=msar[:,cropseq[0]:cropseq[1]]
    
    sl=len(msar[0].seq)

    fn=shade.seqfeat2shadefeat(msar,feature_types=feature_types,force_feature_pos='bottom',debug=debug)
    shaded=ipyshade.shadedmsa4plot(msar,features=fn,shading_modes=shading_modes,debug=debug,startnumber=startnumber,setends=[startnumber-2,sl+startnumber+2])
    
#     shaded=ipyshade.shadepdbquick(pdb_chain_id=pdb_chain_id,entrez_email='info@example.com',debug=False,\
#                          force_feature_pos='bottom',ruler='bottom',legend=False,\
#                          feature_types=['SecStr'],show_seq_names=False,show_seq_length=False,density=dpi,startnumber=1)
    
    #If sl%10=10 se will have a ruler number hanging beyond the sequence image, and we need to correct for that.
    if right_overhang_fix is None:
        if sl%10==0:
            if sl<100:
                rof= 0.1
            else:
                rof=0.5
        else:
            rof=0
    else:
        rof=right_overhang_fix
        
    plot=(ggplot(data=datafixed,mapping=aes(x='resid', y=column))
        + geom_point(size=0.1)+geom_bar(stat='identity',width=0.5)
        + scale_x_continuous(limits=(0.5,sl+0.5+rof),expand=(0,0.2),name='',breaks=[])
       # + scale_y_continuous(breaks=[0,0.5,1.0])
        + theme_light()+theme(aspect_ratio=0.15,dpi=dpi,plot_margin=0)) #+ facet_wrap('~ segid',dir='v')

       
    plot = plot + geom_seq_x(seqimg=shaded.img,xlim=(1,sl+rof),\
                             ylim=(ymin,data[column].max()),aspect_ratio=0.15)
    
    
    return plot
        
        
        


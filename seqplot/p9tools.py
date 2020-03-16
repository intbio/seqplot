from plotnine.utils import resolution
from plotnine.doctools import document
from plotnine.geoms.geom_tile import geom_tile
from plotnine import aes, theme

import pandas as pd
import numpy as np

import matplotlib.image as mpimg
import io
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

@document
class geom_seq_x(geom_tile):
    """
    Rectangles specified using a center points

    {usage}

    Parameters
    ----------
    {common_parameters}

    See Also
    --------
    plotnine.geoms.geom_rect
    """
    DEFAULT_AES = {'alpha': 1, 'color': None, 'fill': '#333333',
                   'linetype': 'solid', 'size': 0.1}
    REQUIRED_AES = {'x', 'y'}
    DEFAULT_PARAMS = {'stat': 'identity', 'position': 'identity',
                      'na_rm': False}

    
    def __init__(self, seqimg=None, x=None, y=None, aspect_ratio=0.1, ymin=None, mapping=None, data=None, **kwargs):

        d=mpimg.imread(io.BytesIO(seqimg))
        pict=pd.DataFrame([(i[1],-i[0],to_hex(d[i])) for i in np.ndindex(d.shape[0],d.shape[1])],columns=['x','y','c'])
        imgdx=float(d.shape[1])
        imgdy=float(d.shape[0])
        dy=(y.max()-np.min(y.min(),ymin))
#         dy=103.
        dx=(x.max()-x.min())+1.0
        f=imgdy/imgdx
        a=aspect_ratio
        dh=dy*f/np.abs(a-f)
        pixel_width=dx/imgdx
#         pict.x=pict.x/(imgdx-1)*(dx-pixel_width)-0.5+x.min()+pixel_width/2
        pict.x=pict.x/(imgdx-1)*(dx)-0.5+x.min()
        pict.y=pict.y/imgdy*dh+min(y.min(),ymin)
#         pict.x=pict.x/1416.*103.-0.5
#         pict.y=pict.y/64.*1.

#         super().__init__(aes(x='x', y='yy'),**kwargs)
#         print(self.data)
        super().__init__(aes(x='x', y='y'), pict, fill=pict['c'],**kwargs)

    
    def setup_data(self, data):
        
        
        try:
            width = data.pop('width')
        except KeyError:
            width = resolution(data['x'], False)

        try:
            height = data.pop('height')
        except KeyError:
            height = resolution(data['y'], False)

        data['xmin'] = data['x'] - width / 2
        data['xmax'] = data['x'] + width / 2
        data['ymin'] = data['y'] - height / 2
        data['ymax'] = data['y'] + height / 2
        return data

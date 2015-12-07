"""
    theme.py

    Author: H.C. v. Stockhausen < hc at vst.io >
    Date:   2012-10-14
    
"""

from __future__ import absolute_import
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib import colors, enums
from nougat.pdf.common import *


class DefaultTheme(object):

    _s = getSampleStyleSheet()
    _s2 = getSampleStyleSheet()
    doc = {
        'leftMargin': None,
        'rightMargin': None,
        'topMargin': None,
        'bottomMargin': None
        }
    headers = {
         H1: _s['Heading1'],
         H2: _s['Heading2'],
         H3: _s['Heading3'],
         H4: _s['Heading4'],
         H5: _s['Heading5'],
         H6: _s['Heading6'],
         }
    paragraph = _s['Normal'] 
    paragraph_centered = _s2['Normal']
    paragraph_centered.alignment = enums.TA_CENTER
    spacer_height = 0.25 * inch
    table_style = [
        ('ALIGN', (0,0), (-1,-1), 'LEFT'),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('FONT', (0,0), (-1,0), 'Helvetica-Bold'),
        ('LINEBELOW', (0,0), (-1,0), 1, colors.black),
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#C0C0C0')),
        ('ROWBACKGROUNDS', (0,1), (-1, -1), [colors.white,
        colors.HexColor('#E0E0E0')])]
    
    @classmethod
    def doc_template_args(cls):
        return dict([(k, v) for k, v in cls.doc.items() if v is not None])
    
    @classmethod    
    def header_for_level(cls, level):
        return cls.headers[level]
        
    def __new__(cls, *args, **kwargs):
        raise TypeError("Theme classes may not be instantiated.")
 
 
 

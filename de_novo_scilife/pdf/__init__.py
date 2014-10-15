"""
    __init__.py "

    Author: H.C. v. Stockhausen < hc at vst.io >
    Date:    2012-10-14

"""

from cStringIO import StringIO
import urllib
from reportlab.platypus.doctemplate import SimpleDocTemplate
from reportlab.platypus.flowables import Image
from reportlab.platypus import Paragraph, Spacer, KeepTogether, PageBreak
from reportlab.lib import colors
from reportlab.platypus.tables import Table, TableStyle
from reportlab.platypus import ListFlowable, ListItem
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.rl_config import defaultPageSize
from reportlab.pdfgen import canvas
from theme import DefaultTheme
from util import calc_table_col_widths
from common import *

# Adapted from http://code.activestate.com/recipes/576832/
class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        canvas.Canvas.__init__(self, *args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        """add page info to each page (page x of y)"""
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_number(num_pages)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        self.setFont("Helvetica", 7)
        filename = self._filename.split("/")[-1]
        self.drawCentredString(110*mm, 15*mm,
            "Page %d of %d -- %s" % (self._pageNumber, page_count, filename))


class Pdf(object):

    story = []    
    theme = DefaultTheme
    
    def __init__(self, title, author):
        self.title = title
        self.author = author
    
    def set_theme(self, theme):
        self.theme = theme
        
    def add(self, flowable):
        self.story.append(flowable)
    
    def add_header(self, text, level=H1):
        p = Paragraph(text, self.theme.header_for_level(level))
        self.add(p)
    
    def add_spacer(self, height_inch=None):
        height_inch = height_inch or self.theme.spacer_height
        self.add(Spacer(1, height_inch)) # magic 1? no, first param not yet implemented by rLab guys
        
    def add_paragraph(self, text, style=None):
        style = style or self.theme.paragraph
        p = Paragraph(text, style)
        self.add(p)
    
    def add_list(self, items, list_style=UL):
        styles = getSampleStyleSheet()
        style = styles["Normal"]
        list_to_print = []
        for item in items:
            list_to_print.append(Paragraph(item, style))
        t = ListFlowable(list_to_print , bulletType='i')
        self.add(t)
    
    def add_table(self, rows, width=None, col_widths=None, align=CENTER,
            extra_style=[]):
        style = self.theme.table_style + extra_style
        if width and col_widths is None: # one cannot spec table width in rLab only col widths
            col_widths = calc_table_col_widths(rows, width) # this helper calcs it for us
        table = Table(rows, col_widths, style=style, hAlign=align)
        self.add(table) 
    
    def add_image(self, src, width, height, align=CENTER, caption=None):

        if src.split(".")[-1] in ["png", "PNG"]:
            try:
                f = open(src)
                data = StringIO(f.read())
            except:
                return
            else:
                img = Image(data, width, height)
            finally:
                f.close()
        else:
            img = Image(src, width, height)
        
        img.hAlign = align
        if caption:
            caption_p = Paragraph(caption, self.theme.paragraph_centered)
            image_table = Table([[img], [caption_p]], width)
            image_table.setStyle(TableStyle([('ALIGN',(-1,-1),(-1,-1),'CENTER')]))
            self.add(image_table)
        else:       
            self.add(img)
 
    def add_qrcode(self, data, size=150, align=CENTER):
        "FIXME: ReportLib also supports QR-Codes. Check it out."
        
        src = "http://chart.googleapis.com/chart?"
        src += "chs=%sx%s&" % (size, size)
        src += "cht=qr&"
        src += "chl=" + urllib.quote(data)
        self.add_image(src, size, size, align)

    def add_pagebreak(self):
        self.add(PageBreak())

    def render(self, pdfTitle):
        doc_template_args = self.theme.doc_template_args()
        doc = SimpleDocTemplate("{}".format(pdfTitle), title=self.title, author=self.author,
            **doc_template_args)
        doc.build(self.story, canvasmaker=NumberedCanvas)
    

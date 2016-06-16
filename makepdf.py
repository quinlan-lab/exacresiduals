from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import cm, mm, inch, pica
from reportlab.platypus import Frame, Image
from reportlab.lib import utils

import sys, re

def get_image(path, width=1*cm):
    img = utils.ImageReader(path)
    iw, ih = img.getSize()
    aspect = ih / float(iw)
    return Image(path, width=width, height=(width * aspect))

pdf = Canvas('toppics.pdf', pagesize = letter)
filename=open(sys.argv[1],'r')
story = []
for line in filename:
    frame = Frame(0*inch, 0*inch, 8.5*inch, 11*inch, showBoundary=0)
    fields=re.split('[:-]',line.strip())
    picname='-'.join(fields)
    story.append(get_image('figs/'+picname+'.png', width=8.3*inch))
    frame.addFromList(story, pdf)
    #pdf.drawImage('figs/'+picname+'.png',-.5*inch,0)
    pdf.showPage()
    story = []
pdf.save()

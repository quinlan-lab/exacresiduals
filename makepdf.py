from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import cm, mm, inch, pica
from reportlab.platypus import Frame, Image
from reportlab.lib import utils

import sys, re

from pyPdf import PdfFileReader, PdfFileWriter

filename=open(sys.argv[1],'r')
merged='toppics.pdf'
def create_merged_pdf():
    output = PdfFileWriter()
    for line in filename:
        fields=re.split('[:-]',line.strip())
        picname='-'.join(fields)
        review_pdf = PdfFileReader(file('figs/'+picname+".pdf", "rb"))
        output.addPage(review_pdf.getPage(0))
    outputStream = file(merged, "wb")
    output.write(outputStream)
    outputStream.close()

create_merged_pdf()
 
# def get_image(path, width=1*cm):
    # img = utils.ImageReader(path)
    # iw, ih = img.getSize()
    # aspect = ih / float(iw)
    # return Image(path, width=width, height=(width * aspect))

# pdf = Canvas('toppics.pdf', pagesize = letter)
# filename=open(sys.argv[1],'r')
# story = []
# for line in filename:
    # frame = Frame(0*inch, 0*inch, 8.5*inch, 11*inch, showBoundary=0)
    # fields=re.split('[:-]',line.strip())
    # picname='-'.join(fields)
    # story.append(get_image('figs/'+picname+'.pdf', width=8.3*inch))
    # frame.addFromList(story, pdf)
    # #pdf.drawImage('figs/'+picname+'.png',-.5*inch,0)
    # pdf.showPage()
    # story = []
# pdf.save()

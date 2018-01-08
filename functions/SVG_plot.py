import cairosvg

class SVG_plot:
    '''
    This class contains all the methods and data to create the genome plot in svg format.
    Once the process is done, use can choose to save the svg or render using cairosvg.
    Good luck boy.
    '''

    def __init__(self, width, height, pixel):
        self.width = width
        self.height = height
        self.pixel = pixel
        self.plot =  '<svg width="%s" height="%s">\n' % (width*pixel, height*pixel)

    # Adding square:
    def draw_chunk(self, row):
        '''
        This function expects a row of a dataframe in which there must be the
        following three columns:
            x - column number
            y - row number
            color - color of the field in hexadecimal code.
        '''
        x = row['x']
        y = row['y']
        color = row['color']
        self.plot += ('<rect x="%s" y="%s" width="%s" height="%s" style="stroke-width:1;stroke:%s; fill: %s" />\n' %
            (x*self.pixel, y*self.pixel, self.pixel , self.pixel, color, color))

    # Adding dot:
    def draw_GWAS(self, row):
        x = row['x']
        y = row['y']
        self.plot += ('<circle cx="%s" cy="%s" r="2" stroke="%s" stroke-width="1" fill="%s" />\n' %
            (x*self.pixel, y*self.pixel, "black", "black"))

    def mark_centromere(self, centr_start, centr_end):
        # Calculate centromere dimensions:
        centr = (
            int((centr_start / self.width) * self.pixel),
            int((centr_end / self.width) * self.pixel),
            int((centr_end + centr_start) / (self.width * 2) * self.pixel),
            int((centr_end - centr_start) / (self.width * 2) * self.pixel)
        )

        # Marking centromere on the left:
        self.plot += ('<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>' %
                 ( 0, centr[0],
                   0, centr[2],
                   centr[3], centr[2],
                   centr[3]*2, centr[2],
                   centr[3], centr[2],
                   0, centr[2],
                   0, centr[1]))

        # Marking centromoere on the right:
        self.plot += ('<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>' %
                 ( self.width * self.pixel, centr[0],
                   self.width * self.pixel, centr[2],
                   self.width * self.pixel - centr[3], centr[2],
                   self.width * self.pixel - centr[3]*2, centr[2],
                   self.width * self.pixel - centr[3], centr[2],
                   self.width * self.pixel, centr[2],
                   self.width * self.pixel, centr[1]))

    # Adding legend?
    # Adding centromere?
    # Adding frame?
    # Adding custom annotation.


    # Close svg document:
    def __close(self):
        self.plot += '</svg>\n'

    # Save svg file:
    def save_svg(self, file_name):
        self.__close() # Closing the svg tag.
        f = open(file_name, 'w')
        f.write(self.plot)
        f.close()

    # Render png file:
    def save_png(self, file_name):
        cairosvg.svg2png(bytestring=self.plot,write_to=file_name)

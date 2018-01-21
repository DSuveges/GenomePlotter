# The svg class
import cairosvg
import math
import sys

class SVG_plot:
    '''
    This class contains all the methods and data to create the genome plot in svg format.
    Once the process is done, use can choose to save the svg or render using cairosvg.

    margins = array of four integers corresponding to the left, upper, right and bottom margin.

    At this point the function won't check if the arguments are correct.
    '''

    def __init__(self, width, height, pixel, margins = [0,0,0,0]):

        # Testing if margins are proper:
        try:
            self.margins = [ int(x * pixel) for x in margins]
        except:
            sys.exit("[Error] Margins were not properly formatted! Expecting an array of four floats.")

        self.width = width
        self.height = height
        self.pixel = pixel
        self.margins = margins

        # Constants for the cytobnand plot:
        self.cband_x1 = self.margins[0] * 0.5
        self.cband_x2 = self.margins[0] * 0.7
        self.cband_width = (self.cband_x2 - self.cband_x1)
        self.cband_colors = {
            'gneg'    : '#FFFFFF', # White
            'gpos25'  : '#E5E5E5', # Gray90
            'gpos50'  : '#CCCCCC', # Gray80
            'gpos75'  : '#B3B3B3', # Gray70
            'gpos100' : '#999999', # Gray60
            'acen'    : '#CCCCCC', # Gray80
            'gvar'    : '#999999', # Gray60
            'stalk'   : '#E5E5E5', # Gray90
            'border'  : '#999999'  # Gray60
        }

        # This is a scaled star for plotting custom annotation:
        self.star = [e for ts in zip([self.__rotate(0, pixel*2, math.radians(alpha*36)) for alpha in range(0,9,2)],
                                    [self.__rotate(0, pixel*2, math.radians(alpha*36), 1.7) for alpha in range(1,10,2)],) for e in ts]

        self.plot =  '<svg width="%s" height="%s">\n' % (width*pixel + margins[0] + margins[2], height*pixel + margins[1] + margins[3])

    def __adjust_coord(self, x, y):
        ''' This simple function just sifts the coordinates by the margins and corrects for the pixel size.'''
        return (int(x*self.pixel + self.margins[0]),
                int(y*self.pixel + self.margins[1]))

    def add_assoc(self, row):
        ''' Adding a star to the polt and the rsID and trait... once later'''
        (x,y) = self.__adjust_coord(row['x'], row['y'])

        # Offsetting the star:
        self.plot += ('<polygon stroke="white" fill="#E25FDA" stroke-width="2" points="%s" />\n' %
                      " ".join([",".join([str(x + a[0] + self.pixel/2),str(y+a[1] + self.pixel/2)]) for a in self.star]))

        # Adding annotation:
        # to be implemented.... not yet there...


    # Adding square:
    def draw_chunk(self, row, y_scale = 1, y_shift = 0):
        '''
        This function expects a row of a dataframe in which there must be the
        following three columns:
            x - column number
            y - row number
            color - color of the field in hexadecimal code.

        y_scale: the pre-defined pixel size will be scaled by this factor.
        y_sihft: the box will be shifted on the y axis by this factor (in pixel)
        '''
        (x,y) = self.__adjust_coord(row['x'],row['y'])
        color = row['color']
        self.plot += ('<rect x="%s" y="%s" width="%s" height="%s" style="stroke-width:1;stroke:%s; fill: %s" />\n' %
            (x, (y + self.pixel * y_shift), self.pixel , self.pixel * y_scale, color, color))

    # Adding dot:
    def draw_GWAS(self, row):
        (x,y) = self.__adjust_coord(row['x'],row['y'])
        self.plot += ('<circle cx="%s" cy="%s" r="%s" stroke="%s" stroke-width="1" fill="%s" />\n' %
            (x + self.pixel/2, y + self.pixel/2, self.pixel * 0.75, "black", "black"))

    def mark_centromere(self, centr_start, centr_end):

        # calculating the y coordinates of the centromeres based on the chunk count:
        start_y = int(centr_start / self.width)
        end_y = int(centr_end / self.width)

        # Marking centromere on the left:
        self.plot += ('<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>\n' %
                  ( self.__adjust_coord(0, start_y) +
                    self.__adjust_coord(0, (start_y+end_y)/2) +
                    self.__adjust_coord((end_y-start_y), (start_y+end_y)/2) +
                    self.__adjust_coord((end_y-start_y)*2, (start_y+end_y)/2) +
                    self.__adjust_coord((end_y-start_y), (start_y+end_y)/2) +
                    self.__adjust_coord(0, (start_y+end_y)/2) +
                    self.__adjust_coord(0, end_y )))

        # Marking centromoere on the right:
        self.plot += ('<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>\n' %
                  ( self.__adjust_coord(self.width+1, start_y) +
                    self.__adjust_coord(self.width+1, (start_y+end_y)/2) +
                    self.__adjust_coord(self.width+1-(end_y-start_y), (start_y+end_y)/2) +
                    self.__adjust_coord(self.width+1-(end_y-start_y)*2, (start_y+end_y)/2) +
                    self.__adjust_coord(self.width+1-(end_y-start_y), (start_y+end_y)/2) +
                    self.__adjust_coord(self.width+1, (start_y+end_y)/2) +
                    self.__adjust_coord(self.width+1, end_y)))
    # Adding legend?
    # Adding centromere? DONE
    # Adding frame?
    # Adding custom annotation.
    # Adding cytobands.
    # Rotating and scaling vecot
    def __rotate(self, x, y, angle, scale=1):
        return ((x * math.cos(angle) - y * math.sin(angle))*scale,
                (x * math.sin(angle) + y * math.cos(angle))*scale)

    # Close svg document:
    def __close(self):
        self.plot += '</svg>\n'

    # Drawing the cytoband ruler next to the chromosome:
    def draw_cytoband(self, row):

        (x1,y1) = self.__adjust_coord(row['x'],row['y1'])
        (x2,y2) = self.__adjust_coord(row['x'],row['y2'])

        # Centromeres are plotted as triangles facing against each other:
        if row['type'] == 'acen' and 'p' in row['name']:
            self.plot += ('<polygon points="%s,%s %s,%s %s,%s" style="fill:%s;stroke:%s;stroke-width:1;fill-rule:nonzero;" />\n' %
                  (self.cband_x1, y1, self.cband_x2, y1, (self.cband_x1 + self.cband_x2) / 2 , y2, self.cband_colors[row['type']], self.cband_colors['border']))
        elif row['type'] == 'acen' and 'q' in row['name']:
            self.plot += ('<polygon points="%s,%s %s,%s %s,%s" style="fill:%s;stroke:%s;stroke-width:1;fill-rule:nonzero;" />\n' %
                  (self.cband_x1, y2, self.cband_x2, y2, (self.cband_x1 + self.cband_x2) / 2, y1, self.cband_colors[row['type']], self.cband_colors['border']))

        # Regular bands are plotted as rectables, where the fill color is set based on the type of the band:
        else:
            self.plot += ('<rect x="%s" y="%s" width="%s" height="%s" style="stroke-width:1;stroke:%s; fill: %s" />\n' %
                (self.cband_x1, y1, self.cband_width , y2 - y1, self.cband_colors['border'], self.cband_colors[row['type']]))

        # Adding the name of the band:
        self.plot += ('<text x="%s" y="%s" text-anchor="end" font-family="sans-serif" font-size="50px" fill="%s">%s</text>' %(self.cband_x1*0.8, y1 + (y2 - y1)/2, self.cband_colors['border'], row['name']))

    # Save svg file:
    def save_svg(self, file_name):
        self.__close() # Closing the svg tag.
        f = open(file_name, 'w')
        f.write(self.plot)
        f.close()

    # Render png file:
    def save_png(self, file_name):
        cairosvg.svg2png(bytestring=self.plot,write_to=file_name)

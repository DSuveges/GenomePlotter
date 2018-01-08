
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


    # Adding other?


    # Save svg:


    # Render png:

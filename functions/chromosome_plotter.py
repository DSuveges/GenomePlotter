import pickle
import pandas as pd
import cairosvg



class chromosome_plotter(object):

    def __init__(self, input_data, pixel):
        self.__pixel__ = pixel
        self.__chromosomeData__ = input_data
        self.__chromosomeName__ = input_data.chr[0]
        self.__chunk_size__ = input_data.iloc[0].end

        # Calculate width and height:
        self.__widthChunk__ = input_data.x.max()
        self.__width__ = pixel * (input_data.x.max() + 1)
        self.__height__ = pixel * (input_data.y.max() + 1)

        # Do something:
        self.__plotString__ = ''

    def draw_dummy(self, color = '#B3F29D'):
        self.__plotString__  = ('<rect x="0" y="0" width="%s" height="%s" style="stroke-width:1;stroke:%s; fill: %s" />\n' %
            (self.__width__, self.__height__, color, color))

    def add_centromere(self, cytobandDf):
        # Get start and end positions of the centromere from the cytoband dataframe
        chromosome = self.__chromosomeName__
        centromereStarts = cytobandDf.loc[(cytobandDf.chr == str(chromosome)) & (cytobandDf.type == 'acen') , 'start'].tolist()
        centromereEnds = cytobandDf.loc[(cytobandDf.chr == str(chromosome)) & (cytobandDf.type == 'acen') , 'end'].tolist()
        
        cetromereXCoord = self.__width__ / 3

        # Calculating important measures from the cytobands dataframe:
        centromereStartCoord = int((centromereStarts[0] * self.__pixel__) / (self.__chunk_size__ * self.__widthChunk__))
        centromereEndCoord = int((centromereEnds[1] * self.__pixel__) / (self.__chunk_size__ * self.__widthChunk__))
        centromereMidpoint = (centromereEndCoord + centromereStartCoord) / 2
        centromereMidpointDiff = (centromereStartCoord - centromereStartCoord) / 2 

        # Right mark of the centromere:
        half_centromere = ('<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>\n' %
            (0, centromereStartCoord,  
             0, centromereMidpoint, cetromereXCoord/2, centromereMidpoint, cetromereXCoord, centromereMidpoint, 
             cetromereXCoord/2, centromereMidpoint, 0, centromereMidpoint, 0,centromereEndCoord))

        # Generating the other half of the centromoere:
        other_half = '<g transform="rotate(180 0 %s) translate(-%s, -%s)">\n\t%s\n</g>\n' %(
            centromereMidpoint, self.__width__, 0, half_centromere)

        # adding both sides of the centromere to the plot:
        self.__plotString__ += '<g id="centromere">\n%s %s </g>\n' % (half_centromere, other_half)

    def wrap_svg(self, fileName):
        self.__svg__ =  '<svg width="%s" height="%s">\n%s</svg>' % (self.__width__, self.__height__, self.__plotString__)

        f = open(fileName, 'w')
        f.write(self.__svg__)
        f.close()

    def return_data(self, dataType):
        if dataType == 'svg':
            return(self.__svg__)
        elif dataType == 'data':
            return(self.__chromosomeData__)

    def pickle_data(self, fileName):
        pickle.dump( {
            'chromosome' : self.__chromosomeName__, 
            'width'      : self.__width__,
            'height'     : self.__height__,
            'pixel'      : self.__pixel__,
            'chunk_size' : self.__chunk_size__,
            'data'       : self.__plotString__,
        }, open( fileName, "wb" ) )

    def save_png(self, fileName):
        cairosvg.svg2png(bytestring=self.__svg__,write_to=fileName)

    def draw_chromosome(self):
        self.__plotString__  = ''
        self.__chromosomeData__.head()

        def add_chunk(df_row, pixel):
            color = df_row['color']
            x = df_row['x'] * pixel
            y = df_row['y'] * pixel
            self.__plotString__  += ('<rect x="%s" y="%s" width="%s" height="%s" style="stroke-width:1;stroke:%s; fill: %s" />\n' %
               (x, y, pixel, pixel, color, color))

        self.__chromosomeData__.apply(add_chunk, axis = 1, args = (self.__pixel__,))
        





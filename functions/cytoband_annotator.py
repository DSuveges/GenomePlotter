import pandas as pd
import cairosvg

def get_centromere_position(cytobandFile, chromosome):
    """
    This function parses the gzipped cytoband file and returns the location of the centromere for a given chromosome.
    
    Input:
    cytobandFile = gzipped cytoband data
    chromosome = chromosome name as string.
    
    Output:
    centromere position
    """
    
    # Reading cytoband file as a pandas dataframe:
    df = pd.read_csv(cytobandFile, sep = "\t", compression='gzip')
    
    # Extracting centromere for that chromosome:
    cytobands = df.loc[df.chr == chromosome]
    if len(cytobands) == 0:
        print('[Error] Cytobands were not found for chromosome {}.'.format(chromosome))
        return None
        
    # Extracting centromere:
    centromere = cytobands.loc[cytobands.type == 'acen']
    if len(centromere) == 0:
        print('[Error] Centromeres were not found for chromosome {}.'.format(chromosome))
        return None
    
    # Extracting midpoint:
    centr = centromere.loc[centromere.name.str.match('q')].end.tolist()[0]
    return centr

class cytoband_annotator(object):
    # Built in strings with the cytoband definitions:
    centromer = '<polygon points="{x1},{y1} {x2},{y2} {x3},{y3}" style="fill:{fill_color};stroke:{border_color};stroke-width:{box_width};fill-rule:nonzero;" />\n'
    cytoband_name = '<text x="{x}" y="{y}" text-anchor="end" font-family="sans-serif" font-size="{font_size}px" fill="{font_color}">{band_name}</text>\n'
    cytoband_box = '<rect x="{x}" y="{y}" width="{width}" height="{height}" style="stroke-width:{box_width};stroke:{border_color};fill:{fill_color}" />\n'

    def __init__(self, pixel, chromosome, bandFile, chunkSize, width, cytbandColors):
        '''
        This class generates gwas signals based on the provided parameters
        '''
        
        # Reading gwas file:
        cytobandDf = pd.read_csv(bandFile, compression='gzip', sep='\t', quotechar='"', header=0,)

        # Filtering cytoband dataframe;
        cytobandDf_select = cytobandDf.loc[ cytobandDf.chr == chromosome]

        # Get chunk count from position:
        cytobandDf_select = cytobandDf_select.assign(startY = cytobandDf_select.start / chunkSize / width * pixel)
        self.cytobandDf_select = cytobandDf_select.assign(endY = cytobandDf_select.end / chunkSize/ width * pixel)
        
        # cytoband drawing parameters derived from the plot parameters:
        self.x_offset = pixel * 40
        self.font_size = pixel * 9
        self.border_width = int(pixel / 2)
        self.box_width = pixel * 5
        self.right_margin = pixel * 4
        
        # Color definitions:
        self.cytbandColors = cytbandColors
    
    def generate_bands(self):
        # Shorhands for box parameters:
        x0 = self.x_offset # cytoband box offset
        bw = self.border_width # cytoband border width
        bxw = self.box_width # cytoband box width
        chw = bxw/2+x0 # centromer half point
        
        def _temp(row):
            svg_element = ''

            # Centromeres are plotted as triangles facing against each other:
            if row['type'] == 'acen' and 'p' in row['name']:
                
                svg_element = self.centromer.format(** {'x1' : x0, 'y1' : row['startY'],
                                                  'x2' : x0 + bxw, 'y2' : row['startY'], 
                                                  'x3' : chw, 'y3' : row['endY'], 
                                                  'fill_color' : self.cytbandColors['acen'], 
                                                  'border_color' : self.cytbandColors['border'], 
                                                  'box_width' : bw})
            elif row['type'] == 'acen' and 'q' in row['name']:
                
                svg_element = self.centromer.format(** {'x1' : x0, 'y1' : row['endY'],
                                                  'x2' : x0 + bxw, 'y2' : row['endY'], 
                                                  'x3' : chw, 'y3' : row['startY'], 
                                                  'fill_color' : self.cytbandColors['acen'], 
                                                  'border_color' : self.cytbandColors['border'], 
                                                  'box_width' : bw})
                                                  
            # Regular bands are plotted as rectangles, where the fill color is set based on the type of the band:
            else:
                svg_element = self.cytoband_box.format(** {'x' : x0, 'y' : row['startY'],
                                                    'width' : bxw,
                                                    'height' : row['endY'] - row['startY'],
                                                    'fill_color' : self.cytbandColors[row['type']], 
                                                    'border_color' : self.cytbandColors['border'], 
                                                    'box_width' : bw})

            # Adding cytoband names, except centromeres:
            if row['type'] != 'acen':
                svg_element += self.cytoband_name.format(**{ 'x': x0 * 0.8, 
                                                    'y' : row['startY'] + (row['endY'] - row['startY'])/2,
                                                    'font_size' : self.font_size,
                                                    'font_color' : self.cytbandColors['border'], 
                                                    'band_name' : row['name']})

            return(svg_element)

        self.bands = self.cytobandDf_select.apply(_temp, axis = 1)
        
    def generate_png(self, filename = 'test_box.png'):
        (width, height) = self.get_dimensions()
        
        bands_string = '<svg width="{}" height="{}">'.format(width, height)
        bands_string += '\n'.join(self.bands.tolist())
        bands_string += '</svg>'
        cairosvg.svg2png(bytestring=bands_string,write_to=filename)
        
    def return_svg(self):
        return( '\n'.join(self.bands.tolist()))
    
    def get_dimensions(self):
        width = self.x_offset + self.box_width + self.right_margin
        height = int(self.cytobandDf_select.endY.max())
        return(width, height)


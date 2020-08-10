import pandas as pd


# This helper class converts genomic position to screen position:
class position_converter(object):
    def __init__(self, width, chunkSize, pixel):
        self.width = width
        self.chunkSize = chunkSize
        self.pixel = pixel
        
    def convert(self, position):
        return(position / (self.width * self.chunkSize) * self.pixel)

class gene_annotator(object):
    
    # Text line:
    text = '<text x="{x}" y="{y}" alignment-baseline="baseline" text-anchor="start" font-family="sans-serif" font-size="{font_size}px" fill="{font_color}">{gene_name}</text>\n'
    
    # Segment line:
    segment = '<line x1="{x1}" x2="{x2}" y1="{y1}" y2="{y2}" stroke="black" stroke-width="2"/>\n'
    
    # Font color:
    __fontColor = 'black'
    
    def __init__(self, gene_file, chromosome_svg_data):
        """
        The gene file is supposed to be a compressed tab limited file.
        Expected columns: chr, start, end, name
        """
        
        # Open gene file and store data:
        self.__df = pd.read_csv(gene_file, sep= "\t", compression='gzip')
        
        # Store other values:
        self.__chunkSize = chromosome_svg_data['chunk_size']
        self.__pixel = chromosome_svg_data['pixel']
        self.__svg_top = 0
        self.__svg_bottom = chromosome_svg_data['height']
        self.__svg_width = chromosome_svg_data['width']
        self.__width = chromosome_svg_data['width'] / chromosome_svg_data['pixel']
        self.__fontSize = chromosome_svg_data['pixel'] * 10

        # Initialize return values:
        self.__gene_annot = ''

    def genes_on_chromosome_arms(self, genes_df, centromerePosition, arm):

        # Initialize position converter:
        converter = position_converter(self.__width, self.__chunkSize, self.__pixel)

        # Get centromere position:
        limit = converter.convert(centromerePosition)

        # Extract values:
        font_size = self.__fontSize
        font_color = self.__fontColor
        
        for i, row in genes_df.iterrows():
            pos = converter.convert(row['start'])
            
            if arm == 'p':
                if limit <= pos:
                    text_pos = limit - 10
                    limit = text_pos - font_size - 10
                else:
                    text_pos = pos
                    limit = text_pos - font_size - 10
            if arm == 'q':
                if limit > pos:
                    text_pos = limit + 10
                    limit = text_pos + font_size + 10
                else:
                    text_pos = pos
                    limit = text_pos + font_size + 10
            
            # Update svg_top and bottom if required:
            if limit - 50 < self.__svg_top: 
                self.__svg_top = limit - font_size - 50 

            if limit + 50 > self.__svg_bottom:
                self.__svg_bottom = limit + font_size + 50
                            
            # Add segments:
            self.__gene_annot += self.segment.format(**{
                'x1' : 0,
                'x2' : 100,
                'y1' : pos,
                'y2' : pos
            })
            self.__gene_annot += self.segment.format(**{
                'x1' : 100,
                'x2' : 200,
                'y1' : pos,
                'y2' : text_pos
            })
            self.__gene_annot += self.segment.format(**{
                'x1' : 200,
                'x2' : 300,
                'y1' : text_pos,
                'y2' : text_pos
            })

            # Add text:
            self.__gene_annot += self.text.format(**{ 'x': 200, 
                'y' : text_pos - 10,
                'font_size' : font_size,
                'font_color' : font_color, 
                'gene_name' : row['name']}
            )
                      
    def generate_gene_annotation(self, chromosome, centromerePosition):

        # Initialize gene annotation:
        self.__gene_annot = ''
        
        # Restrict data to the requested chromosome:
        print("[Info] Filtering gene data... ")
        
        genes_df = self.__df.loc[ self.__df.chr == chromosome]
        
        print("[Info] There are {} genes on this chromosome.".format(len(genes_df)))
        
        # Generate text for the p-arm:
        self.genes_on_chromosome_arms(genes_df.loc[ genes_df.start <= centromerePosition].sort_values(by=['start'], ascending=False), centromerePosition, 'p')

        # Generate text for the q-arm:
        self.genes_on_chromosome_arms(genes_df.loc[ genes_df.start > centromerePosition].sort_values(by=['start'], ascending=True), centromerePosition, 'q')


    def get_annotation(self):
        return self.__gene_annot

    def get_dimensions(self):
        return (self.__svg_top, self.__svg_bottom)

import pickle
import pandas as pd
import cairosvg
from tqdm import tqdm


class chromosome_plotter(object):

    chunk_svg = '<rect x="{}" y="{}" width="{}" height="{}" style="stroke-width:1;stroke:{}; fill: {}" />'

    def __init__(self, input_data, pixel):
        self.__pixel__ = pixel
        self.__chromosomeData__ = input_data
        self.__chromosomeName__ = input_data.chr[0]
        self.__chunk_size__ = input_data.iloc[0].end

        # Calculate width and height:
        self.__widthChunk__ = input_data.x.max()
        self.__width__ = pixel * (input_data.x.max() + 1)
        self.__height__ = pixel * (input_data.y.max() + 1)

        # The svg string will be stored here:
        self.__plotString__ = ''

    def __add_centromere(self):

        # We use the Gencode annotation in the chromosome dataframe to get start and end:
        centromere_start = self.__chromosomeData__.loc[ self.__chromosomeData__.GENCODE == 'centromere' ].y.min() * self.__pixel__
        centromere_end = self.__chromosomeData__.loc[ self.__chromosomeData__.GENCODE == 'centromere' ].y.max() * self.__pixel__ 
        
        centromere_midpoint = (centromere_end - centromere_start)/2
        centromere_hight = centromere_end - centromere_start


        # How deep we want the cleavage:
        cetromere_x = self.__width__ / 3

        # Right mark of the centromere:
        centromere_string = (f'<path d="M -1 0 \
            C 0 {centromere_midpoint}, {cetromere_x/2} {centromere_midpoint}, {cetromere_x/2} {centromere_midpoint} \
            C {cetromere_x/2} {centromere_midpoint}, 0 {centromere_midpoint}, -1 {centromere_hight} Z" fill="white"/>\n' )


        half_centromere = f'<g transform="translate(0, {centromere_start})">\n\t{centromere_string}\n</g>\n'

        # Generating the other half of the centromoere:
        other_half = f'<g transform="rotate(180 0 {centromere_start + centromere_midpoint}) translate(-{self.__width__}, 0)">\n\t{half_centromere}\n</g>\n'

        # adding both sides of the centromere to the plot:
        self.__plotString__ += f'\n<g id="centromere">\n\t{half_centromere}\t{other_half}</g>\n'


    def draw_dummy(self):
        width = self.__width__
        height = self.__height__

        # Extract dummy and centromere color:
        dummy_color = self.__chromosomeData__.loc[ self.__chromosomeData__.GENCODE != 'centromere' ].color.tolist()[0]
        centromere_color = self.__chromosomeData__.loc[ self.__chromosomeData__.GENCODE == 'centromere' ].color.tolist()[0]
        
        # Extract centromere positions:
        centromere_start = self.__chromosomeData__.loc[ self.__chromosomeData__.GENCODE == 'centromere' ].y.min() * self.__pixel__
        centromere_end =  self.__chromosomeData__.loc[ self.__chromosomeData__.GENCODE == 'centromere' ].y.max() * self.__pixel__ - centromere_start

        print(f'centromere_start: {centromere_start}, centromere_end: {centromere_end}')

        # Adding the full chromosome in dummy;
        self.__plotString__ += self.chunk_svg.format(0, 0, width, height, dummy_color, dummy_color)

        # Adding centromere rectangle:
        self.__plotString__ += self.chunk_svg.format(0, centromere_start, width, centromere_end, centromere_color, centromere_color)

        # Adding centromoere:
        self.__add_centromere()


    def draw_chromosome(self):
        
        pixel = self.__pixel__
        svg_chunks = []
        for index, df_row in  tqdm(self.__chromosomeData__.iterrows(), total=self.__chromosomeData__.shape[0]):
            x = df_row['x'] * pixel
            y = df_row['y'] * pixel
            svg_chunks.append(self.chunk_svg.format(x, y, pixel, pixel, df_row['color'], df_row['color']))

        self.__plotString__ = '\n'.join(svg_chunks)

        # Adding centromoere:
        self.__add_centromere()


    
    def get_plot_with(self):
        return self.__width__


    def get_plot_height(self):
        return self.__height__


    def return_svg(self):
        return self.__plotString__


    def save_png(self, fileName):
        cairosvg.svg2png(bytestring=self.__svg__,write_to=fileName)


    def wrap_svg(self, fileName):
        self.__svg__ =  '<svg width="%s" height="%s" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve">\n%s</svg>' % (self.__width__, self.__height__, self.__plotString__)

        f = open(fileName, 'w')
        f.write(self.__svg__)
        f.close()



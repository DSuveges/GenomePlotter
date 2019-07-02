import pandas as pd

class gwas_annotator(object):
    def __init__(self, pixel, chromosome, gwasFile, chunkSize, width, xoffset = 0, yoffset = 0, gwasColor = "black"):
        '''
        This class generates gwas signals based on the provided parameters
        '''
        
        # Reading gwas file:
        gwas_df = pd.read_csv(gwasFile, compression='gzip', sep='\t', quotechar='"', header=0, 
                              dtype={'#chr': str, 'start' : int, 'end' : int, 'rsID': str, 'trait' : str})

        # Filtering gwas file:
        filtered_locations = gwas_df.loc[gwas_df['#chr'] == chromosome,'start'].tolist()
        
        # Calculate positions:
        chunks = [ int( x / chunkSize) for x in filtered_locations ]
        chunks = set(chunks)
        
        # Based on the chunk number, we calculate the xy position on the 
        self.__positions = [ (int(x % width), int(x / width) ) for x in chunks]
        
        self.__pixel = pixel
        self.__xoffset = xoffset
        self.__yoffset = yoffset
        self.__gwasColor = gwasColor
    
    def generateGWAS(self):
        pixel = self.__pixel
        positions = self.__positions
        xoffset = self.__xoffset
        yoffset = self.__yoffset
        gwasColor = self.__gwasColor
        
        # GWAS Points:
        gwasPoints = ''
        
        # Based on the x/y coordinates, let's draw the point:
        for position in positions:
            gwasPoints += ('<circle cx="%s" cy="%s" r="%s" stroke="%s" stroke-width="1" fill="%s" />\n' %(
                position[0] * pixel + pixel / 2 + xoffset, 
                position[1] * pixel + pixel / 2 + yoffset, 
                pixel/2,
                gwasColor, gwasColor
            ))

        return(gwasPoints)
    
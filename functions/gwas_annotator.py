import pandas as pd
import math

class gwas_annotator(object):

    """Adds GWAS associations to the chromosome"""

    # GWAS hit svg definition:
    gwas_hit = '<circle cx="{}" cy="{}" r="{}" stroke="{}" stroke-width="1" fill="{}" />'

    # Unit circle:
    circle_unit = 63

    # GWAS hit count cap:
    gwas_cap = 10

    def __init__(self, pixel, chromosome, gwasFile, chunkSize, width, xoffset=0, yoffset=0, gwasColor="black"):
        """
        This class generates gwas signals based on the provided parameters
        """

        # Reading gwas file:
        gwas_df = pd.read_csv(gwasFile, compression='gzip', sep='\t', quotechar='"', header=0,
                              dtype={'#chr': str, 'start': int, 'end': int, 'rsID': str, 'trait': str})

        # Filtering dataframe for the given chromosome:
        filtered_locations = gwas_df.loc[gwas_df['#chr'] == chromosome]

        # Looping through all GWAS hits on the chromosome and calculate coordinates and scale:
        data = []
        for index, value in filtered_locations.start.apply(lambda x: int(x / chunkSize)).value_counts().iteritems():
            value = value if value < 10 else 10
            data.append({
                'counts': value,
                'x': int(index % width),
                'y': int(index / width)
            })

        self.__positions = pd.DataFrame(data)

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
        gwasPoints = []

        # Based on the x/y coordinates, let's draw the point:
        for i, row in positions.iterrows():

            # The radius of the circle is proportional to the number of GWAS hits in the given chunk:
            radius = math.sqrt(row['counts'] * self.circle_unit / math.pi)

            # Adding point{}
            gwasPoints.append(self.gwas_hit.format(
                (row['x'] * pixel) + radius / 2 + xoffset,
                (row['y'] * pixel) + radius / 2 + yoffset,
                radius,
                gwasColor, gwasColor
            ))

        return '\n'.join(gwasPoints)

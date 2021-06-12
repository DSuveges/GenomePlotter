import gzip
class process_chrom:
    """This class plots the chromosomes and saves a gzipped text file."""

    def __init__(self, width, height, pixel, defs=0):
        self.width = width
        self.height = height
        self.pixel = pixel

        # If we would like, we can include the defs tag as well:
        if defs == 1:
            self.__add_def(
                'chunk', '<rect x="0" y="0" width="%s" height="%s" style="stroke-width:0" />\n' %
                (pixel, pixel)
            )

    def __add_def(self, ID, svg_string):
        """
        This function adds an other <g> object to the <defs> section of the
        svg object.

        Until the whole object is saved, all the definitions will be stored in a
        dictionary.
        """

        if 'defs' in self:
            self['defs'][ID] = svg_string
        else:
            self['defs'] = {}
            self['defs'][ID] = svg_string

    # Adding square:
    def draw_chunk(self, row):
        """
        This function expects a row of a dataframe in which there must be the
        following three columns:
            x - column number
            y - row number
            color - color of the field in hexadecimal code.
        """
        x = row['x']
        y = row['y']
        color = row['color']
        self.plot += (
            f'<use x="{x*self.pixel}" y="{y*self.pixel}" href="#chunk" style="fill: {color};"/>\n'
        )

    # Adding dot:
    def draw_GWAS(self, row):
        """
        The GWAS signals are added as dots. There's no need to use g and use tags.
        """
        x = row['x']
        y = row['y']
        self.plot += (
            f'<circle cx="{x*self.pixel}" cy="{y*self.pixel}" r="2" stroke="black" stroke-width="1" fill="black" />\n'
        )

    def mark_centromere(self, centr_start, centr_end):
        """
        Based on the coordinates of the centromoere, this function compiles an svg
        path for the centromere and addes to the plot object.
        """

        # The y coordinate of the centromere start and end points 0 adjusted:
        centr_y_start = 0
        centr_y_end = self.pixel * (centr_end - centr_start) / self.width
        centr_y_midpoint = centr_y_end / 2
        centr_x_midpoint = centr_y_end * 2
        print(centr_y_end, centr_y_start)

        # Creating the centromere path:
        half_centromere = ('<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>' % (
            0, centr_y_start,
            0, centr_y_midpoint,
            centr_x_midpoint / 2, centr_y_midpoint,
            centr_x_midpoint, centr_y_midpoint,
            centr_x_midpoint / 2, centr_y_midpoint,
            0, centr_y_midpoint,
            0, centr_y_end
        ))

        # Set position for the left side:
        left_centromere_side = (
            '<g transform="translate(0, %s)">\n\t%s\n</g>\n' % (
                (centr_start * self.pixel) / self.width, half_centromere
            )
        )

        # Set position for the right side:
        right_centromere_side = '<g transform="rotate(180 0 %s) translate(-%s, -%s)">\n\t%s\n</g>\n' % (
            centr_y_midpoint, self.pixel * self.width, (centr_start * self.pixel) / self.width, half_centromere)

        # Adding centromeres to the plot:
        self.plot += '<g id="centromere">\n%s %s </g>\n' % (left_centromere_side, right_centromere_side)

    # Close svg document:
    def save(self, chr_name):
        """
        This function saves the assembled lines into a textfile.
        We don't expect to be included anything else.w
        """
        with gzip.open(f'chr{chr_name}_genome_chunks.txt.gz', 'wt') as f:
            f.write(self.plot)

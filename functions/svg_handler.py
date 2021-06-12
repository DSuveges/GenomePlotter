import cairosvg

class svg_handler(object):
    """Functions to manipulate svg"""

    __svg_rect__ = '<rect x="{}" y="{}" width="{}" height="{}" style="stroke-width:1;stroke:{}; fill: {}" />\n'
    __svg_label__ = '<text x="{}" y="{}" text-anchor="{}" font-family="sans-serif" \
        font-size="{}px" fill="{}">{}</text>\n'
    __svg_line__ = '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="{}" {} />\n'

    def __init__(self, svg_string, width, height):
        self.__svg__ = svg_string
        self.__width__ = width
        self.__height__ = height

    def group(self, translate=(0, 0)):
        """
        Grouping and transforming object. Should have been more functionally rich
        """
        self.__svg__ = '<g transform="translate(%s %s)">\n%s\n</g>\n' % (
            translate[0], translate[1], self.__svg__)

        # Updating coordinates:
        self.__width__ += abs(translate[0])
        self.__height__ += abs(translate[1])

    def appendSvg(self, svg_string):
        self.__svg__ += svg_string

    def mergeSvg(self, svg_obj):
        self.__svg__ += svg_obj.getSvg()

        self.__width__ = max(self.__width__, svg_obj.getWidth())
        self.__height__ = max(self.__height__, svg_obj.getHeight())

    def __closeSvg(self):
        svg_header = (
            f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" \
            xmlns:xlink="http://www.w3.org/1999/xlink" \
            xml:space="preserve"  width = "{self.__width__}" height = "{self.__height__}">\n'
        )

        self.__closedSVG__ = svg_header + self.__svg__ + '\n</svg>\n'

    def savePng(self, filename='test.png'):
        self.__closeSvg()

        cairosvg.svg2png(bytestring=self.__closedSVG__, write_to=filename)

    def saveSvg(self, filename='test.svg'):
        self.__closeSvg()

        with open(filename, 'w') as f:
            f.write(self.__closedSVG__)

    def getSvg(self):
        return(self.__svg__)

    def getWidth(self):
        return(self.__width__)

    def getHeight(self):
        return(self.__height__)

    def draw_rectangle(self, x, y, width, height, stroke, fill):
        self.__svg__ += self.__svg_rect__.format(x, y, width, height, stroke, fill)

    def draw_line(self, x1, y1, x2, y2, stroke="#000000", stroke_width=3, **kwargs):
        extra_args = ''

        if kwargs:
            for key, value in kwargs.items():
                extra_args += f' {key.replace("_","-")}="{value}"'
        print(extra_args)
        self.__svg__ += self.__svg_line__.format(x1, y1, x2, y2, stroke, stroke_width, extra_args)

    def add_text(self, x, y, text, size=10, fill="#000000", anchor='start'):
        self.__svg__ += self.__svg_label__.format(x, y, anchor, size, fill, text)

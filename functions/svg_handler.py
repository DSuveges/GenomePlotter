import cairosvg

class svg_handler(object):
    def __init__(self, svg_string, width, height):
        self.__svg__ = svg_string
        self.__width__ = width
        self.__height__ = height
        
    def group(self, translate = (0,0)):
        '''
        Grouping and transforming object. Should have been more functionally rich
        '''
        self.__svg__ = '<g transform="translate(%s %s)">\n%s\n</g>\n' % (
            translate[0], translate[1], self.__svg__)
        
        # Updating coordinates:
        self.__width__ += abs(translate[0])
        self.__height__ += abs(translate[1])
        
    def appendSvg(self, svg_string):
        self.__svg__  += svg_string
    
    def mergeSvg(self, svg_obj):
        self.__svg__ += svg_obj.getSvg()
        
        self.__width__ = max(self.__width__, svg_obj.getWidth())
        self.__height__ = max(self.__height__, svg_obj.getHeight())
        
    def __closeSvg(self):
        self.__closedSVG__ = '<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve"  width = "%s" height = "%s">\n%s\n</svg>\n' % (
            self.__width__, self.__height__, self.__svg__)
        
    def savePng(self, filename = 'test.png'):
        self.__closeSvg()

        cairosvg.svg2png(bytestring = self.__closedSVG__, write_to = filename)
    
    def saveSvg(self, filename = 'test.svg'):
        self.__closeSvg()
        
        with open(filename, 'w') as f:
            f.write(self.__closedSVG__)

    def getSvg(self):
        return(self.__svg__)
    
    def getWidth(self):
        return(self.__width__)

    def getHeight(self):
        return(self.__height__)
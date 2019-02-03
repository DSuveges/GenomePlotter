import json 

class readConfig(object):

    def __init__(self, configFileName = "config.json"):
        # Reading file into a document:

        with open(configFileName) as f:
            data = json.load(f)


        self.__data = data;

    def getColorScheme(self):
        return(self.__data['colorScheme'])

    def getFiles(self):
        return(self.__data['inputFiles'])

    def updateConfig(self, configFileNA):
        return(self)


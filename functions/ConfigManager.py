import json 
import os

class ConfigManager(object):
    """
    This class manages all configuration setting: 
    * Set new values
    * Retrives set ones
    * Updates configuration file.
    * Upon retrieveal file and folder configuration, check performed for existence
    """

    def __init__(self, config_file_name):

        # Test if the provided file exists or not:
        if not os.path.isfile(config_file_name):
            raise ValueError(f'The provided config file ({config_file_name}) doesn\'t exits.')

        # Reading json file:
        with open(config_file_name) as f:
            data = json.load(f)

        self.__data = data;
        self.__config_file_name = config_file_name

        # data folder has a special importance:
        self.__data_folder = data['basic_parameters']['data_folder']

        # souces are extracted too:
        self.__sources = data['source_data']

        
    ##
    ## Updating config file
    ##
    def save_config(self, filename=None):

        if not filename:
            filename = self.__config_file_name

        # Saving file:
        with open(filename, 'w') as f:
            f.write(json.dumps(self.__data, indent=4))


    ##
    ## Set functions:
    ##
    def set_data_folder(self, data_folder):
        self.__data_folder = data_folder
        self.__data['basic_parameters']['data_folder'] = data_folder
        

    def set_width(self, width):
        self.__data['plot_parameters']['width'] = width

    
    def set_pixel_size(self, pixel_size):
        self.__data['plot_parameters']['pixel_size'] = pixel_size
        
    
    def set_dark_start(self, dark_start):
        self.__data['plot_parameters']['dark_start'] = dark_start

    
    def set_dark_max(self, dark_max):
        self.__data['plot_parameters']['dark_max'] = dark_max        
    
    
    def set_plot_folder(self, plot_folder):
        self.__data['basic_parameters']['plot_folder'] = plot_folder   


    ##
    ## Get functions:
    ##
    def get_color_scheme(self):
        return self.__data['color_schema']

    
    def get_data_folder(self):
        if not self.__data_folder:
            raise ValueError('Data folder is not yet set!')
        elif not os.path.isdir(self.__data_folder):
            raise ValueError(f"Data folder ({self.__data_folder}) doesn\'t exists.")

        return self.__data_folder

    
    # Functions returning files, tests if the file exists or not:
    def get_source(self, resource, key, chromosome=None):
        """
        This function returns values belonging to a specified key in a selected resource.

        returns: str or int
        """

        if resource not in self.__sources:
            raise ValueError(f'Unknown resource presented: {resource}. Accepted resources: {",".join(self.__sources.keys())}')

        if key not in self.__sources[resource]:
            raise ValueError(f'{key} is not a stored feature for {resource}. Available features: {",".join(self.__sources[resource].keys())}')
            
        value = self.__sources[resource][key]
        
        if chromosome:
            value = value.format(chromosome)
            
        return value
            
        
    def get_cytoband_file(self):
        file = self.__sources['cytoband_data']['processed_file']
        full_path = f'{self.__data_folder}/{file}'
        
        if not os.path.isfile(full_path):
            raise ValueError(f'Cytological band file ({full_path}) doesn\'t exists.')
            
        return full_path
    
    
    def get_chromosome_file(self, chromsome):
        file = self.__sources['ensembl_data']['processed_file']
        file = file.format(chromsome)
        full_path = f'{self.__data_folder}/{file}'
        
        if not os.path.isfile(full_path):
            raise ValueError(f'The requested genome file ({full_path}) doesn\'t exists.')
            
        return full_path    
    
    
    def get_gwas_file(self):
        file = self.__sources['gwas_data']['processed_file']
        full_path = f'{self.__data_folder}/{file}'
        
        if not os.path.isfile(full_path):
            raise ValueError(f'Processed GWAS file ({full_path}) doesn\'t exists.')
            
        return full_path        

    
    def get_gencode_file(self):
        file = self.__sources['gencode_data']['processed_file']
        full_path = f'{self.__data_folder}/{file}'
        
        if not os.path.isfile(full_path):
            raise ValueError(f'Processed GENCODE file ({full_path}) doesn\'t exists.')
            
        return full_path
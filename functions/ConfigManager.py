"""Managing configuration settings for the Genome Plotter."""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from typing import Optional

# List of dataclasses to describe the configuration file:

# List of dataclasses to describe the configuration file:
@dataclass
class PlotParameters:
    """Dataclass to store plot parameters."""
    width: int
    pixel_size: int
    dark_start: float
    dark_max: float
    custom_gene_window: int


@dataclass
class BasicParameters:
    """Dataclass to store basic parameters."""
    chunk_size: int
    missing_tolerance: Optional[float] = None
    plot_folder: Optional[str] = None
    data_folder: Optional[str] = None
    row_length: Optional[int] = None

@dataclass
class ColorSchema:
    """Dataclass to store color schema."""
    chromosome_colors: dict
    cytoband_colors: dict
    gwas_point: str
    arrow_colors: dict


@dataclass
class SourcePrototype:
    """Dataclass to store source prototype."""
    # Output file is mandatory:
    processed_file: str
    # Optional parameters:
    url: Optional[str] = None
    genome_build: Optional[str] = None
    host: Optional[str] = None
    path: Optional[str] = None
    source_file: Optional[str] = None
    arrow_file: Optional[str] = None
    release_date: Optional[str] = None
    release: Optional[int] = None
    version_url: Optional[str] = None
    version: Optional[int] = None

@dataclass
class CytoBandData:
    """Dataclass to store cytoband data."""
    url: str
    processed_file: str
    genome_build: Optional[str] = None

@dataclass
class SourceData:
    """Dataclass to store source data."""
    cytoband_data: CytoBandData
    ensembl_data: SourcePrototype
    gencode_data: SourcePrototype
    gwas_data: SourcePrototype

    def __post_init__(self):
        # Iterating over expected fields and populate data:
        for field in self.__dataclass_fields__.keys():
            self.__setattr__(
                field,
                eval(self.__dataclass_fields__[field].type)(
                    **self.__getattribute__(field)
                )
            )


@dataclass
class Config:
    """Dataclass to store the configuration file."""
    plot_parameters: PlotParameters
    basic_parameters: BasicParameters
    color_schema: ColorSchema
    source_data: SourceData

    def __post_init__(self):
        # Iterating over expected fields and populate data:
        for field in self.__dataclass_fields__.keys():
            class_name =  getattr(
                Config, 
                self.__dataclass_fields__[field].type
            )
            self.__setattr__(
                field,
                class_name(
                    **self.__getattribute__(field)
                )
            )

    # Saving the configuration file:
    def save(self, file_path: str) -> None:
        """Save the configuration file.

        Args:
            file_path (str): Path to the configuration file.
        """
        with open(file_path, "w") as file:
            json.dump(asdict(self), file, indent=3)

    # Updating basic configuration based on command line arguments:
    def update_basic_parameters(self, **kwargs) -> None:
        """Update basic parameters based on command line arguments.

        Args:
            **kwargs: Command line arguments.
        """
        for key, value in kwargs.items():
            if key in self.basic_parameters.__dataclass_fields__.keys():
                self.basic_parameters.__setattr__(key, value)


@dataclass
class ConfigManager(object):
    """
    This class manages all configuration setting:

    * Set new values
    * Retrives set ones
    * Updates configuration file.
    * Upon retrieveal file and folder configuration, check performed for existence
    """

    # Dataclass objects:
    plot_parameters: PlotParameters
    basic_parameters: BasicParameters
    color_schema: ColorSchema
    source_data: SourceData

    def __init__(self, config_file_name):

        # Test if the provided file exists or not:
        if not os.path.isfile(config_file_name):
            raise ValueError(f'The provided config file ({config_file_name}) doesn\'t exits.')

        # Reading json file:
        with open(config_file_name) as f:
            data = json.load(f)

        self.__data = data
        self.__config_file_name = config_file_name

        # data folder has a special importance:
        self.__data_folder = data['basic_parameters']['data_folder']

        # souces are extracted too:
        self.__sources = data['source_data']

    # Updating config file
    def save_config(self, filename=None):

        if not filename:
            filename = self.__config_file_name

        # Saving file:
        with open(filename, 'w') as f:
            f.write(json.dumps(self.__data, indent=4))

    # Set functions:
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

    # Get colors:
    def get_chromosome_colors(self):
        return self.__data['color_schema']['chromosome_colors']

    def get_cytobanc_colors(self):
        return self.__data['color_schema']['cytoband_colors']

    def get_gwas_color(self):
        return self.__data['color_schema']['gwas_point']

    def get_arrow_colors(self):
        return self.__data['color_schema']['arrow_colors']

    def get_color_schema(self):
        return self.__data['color_schema']

    # Get operational parameters:
    def get_data_folder(self):
        if not self.__data_folder:
            raise ValueError('Data folder is not yet set!')
        elif not os.path.isdir(self.__data_folder):
            raise ValueError(f"Data folder ({self.__data_folder}) doesn\'t exists.")

        return self.__data_folder

    def get_source(self, resource, key, chromosome=None):
        """
        This function returns values belonging to a specified key in a selected resource.

        returns: str or int
        """

        if resource not in self.__sources:
            accepted_resources = ",".join(self.__sources.keys())
            raise ValueError(f'Unknown resource presented: {resource}. Accepted resources: {accepted_resources}')

        if key not in self.__sources[resource]:
            available_features = ",".join(self.__sources[resource].keys())
            raise ValueError(f'{key} is not a stored feature for {resource}. Available features: {available_features}')

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

    def get_gencode_arrow_file(self):
        file = self.__sources['gencode_data']['arrow_file']
        full_path = f'{self.__data_folder}/{file}'

        if not os.path.isfile(full_path):
            raise ValueError(f'Processed GENCODE file ({full_path}) doesn\'t exists.')

        return full_path

    # Get plot parameters:
    def get_pixel(self):
        return self.__data['plot_parameters']['pixel_size']

    def get_chunk_size(self):
        return self.__data['basic_parameters']['chunk_size']

    def get_width(self):
        return self.__data['plot_parameters']['width']

    def get_dark_max(self):
        return self.__data['plot_parameters']['dark_max']

    def get_dark_start(self):
        return self.__data['plot_parameters']['dark_start']

    def get_custom_gene_window(self):
        return self.__data['plot_parameters']['custom_gene_window']


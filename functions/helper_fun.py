# Helper functions for the genome plotter:
from time import gmtime, strftime

def generate_xy(df, min_pos, chunk_size, width, position_column='start', x='x', y='y'):

    '''
    Generating x,y coordinates from genomic position and
    the provided plot width or chunk size.

    pos - genomic position
    pixel - the size of the point taken up by one chunk.
    width - number of chunks in one row
    chunk_size - the number of basepairs pooled together in one chunk.
    '''

    pos = int(df[position_column]) - min_pos
    chunk = int(int(pos) / chunk_size)
    df[x] = (chunk - int(chunk / width) * width)
    df[y] = int(chunk / width)
    return df

def get_now():
    return strftime("%H:%M:%S", gmtime())

from fileprocessingsystem import HDFProcessingSystem
from converter import *
import logging
import config

if __name__ == "__main__":

    logging.basicConfig(filename="test.log",filemode="a",format='%(asctime)s %(message)s',level=logging.INFO)
    console = logging.StreamHandler()
    logging.getLogger(config.logger_name).addHandler(console)
    console.setLevel(logging.WARNING)


    file_proc = HDFProcessingSystem()
    C = Converter(input_directory='./data', genetics_filenames = 'FOR2107_22.out.gz', output_filename = 'FOR22b', output_path = './data', chromosome_numbers=22,
                  fileprocessingsystem=file_proc, snp_stats_path = './data/snp-stats')
    output_filepath = os.path.join(C.output_path, C.output_filename + '.h5')
    C.fileprocessingsystem.create_file(output_filepath, 'Genetics') 
    C.run_converter(C.genetics_filenames, ch_nr = '22', hdf_filepath=output_filepath, snp_stats_stem = 'FOR2107_')
    #C.read_snp_stats('./data/snp-stats/FOR2107_22.snp-stats.gz')
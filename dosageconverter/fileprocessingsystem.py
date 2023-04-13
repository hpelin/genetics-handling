from dataclasses import dataclass, field
from converter import FileProcessingSystem
import h5py
import numpy as np


@dataclass
class HDFProcessingSystem(FileProcessingSystem):
    """
    This is a file processing system """
    def read_file(self, file_extension: str, filepath: str):
        return super().read_file(file_extension=file_extension, filepath=filepath)
    

    def save_file(self, hdf_filepath, hdf_group_name, hdf_dataset_name, hdf_map_file_name, dosage_dataset, map_dataset):
        with h5py.File(hdf_filepath, 'a') as hdf:
            G = hdf.create_group(hdf_group_name)
            G.create_dataset(hdf_map_file_name, data=map_dataset, dtype=h5py.string_dtype(), compression = 'lzf')
            G.create_dataset(hdf_dataset_name, data=dosage_dataset, compression = 'lzf')


    def create_file(self, path: str, folder_name: str):
        dt = h5py.special_dtype(vlen=str) 
        with h5py.File(path, 'a') as hdf:
            if folder_name in hdf.keys():
                pass
            else:
                G = hdf.create_group('Genetics')
                G.attrs['map_file_attr'] = np.array(['chromosome', 'snp_id', 'position', 'dosage_ref_all_relative_to', 'refall', 'minor_all_snp_stats', 'major_all_snp_stats','maf','hwe','information'], dtype=dt)

    def check_ch_processed(self, hdf_filepath, hdf_group_name):
        with h5py.File(hdf_filepath, 'a') as hdf:
            if hdf_group_name in hdf.keys():
                return True
            else:
                return False
                           

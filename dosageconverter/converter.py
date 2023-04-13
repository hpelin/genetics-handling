from dataclasses import dataclass, field
import numpy as np
import logging
import os, sys
import gzip 
import utils
from abc import ABC, abstractmethod
import time
import h5py
import config
import string

class FileProcessingSystem(ABC):

    @abstractmethod
    def read_file(self, file_extension: str, filepath: str):
        if file_extension == 'gz':
            read = gzip.open(filepath, 'rt')
            read.seek(0)
        else:
            read = open(filepath, 'r')
            read.seek(0)
        return read
    
    @abstractmethod
    def save_file(self, path: str):
        pass

    def create_file(self, path: str):
        pass


@dataclass
class Converter():
    """
    This is a Converter class containing essentially all the methods 
    needed for the dosage conversion
    """
    fileprocessingsystem: FileProcessingSystem
    showProgress: bool = False
    header: int = 5
    reference_allelle: int = 0
    input_directory: str = ''
    genetics_filenames: list[str] = field(default_factory=list)
    output_path: str = ''
    chromosome_numbers: list[int] = field(default_factory=list)
    dosage_relative_to: int = 3
    process_sample_file: bool = False
    map_file:  list = field(default_factory=list)
    output_filename: str = ''
    logger = logging.getLogger(config.logger_name)
    snp_stats_path: str = ''

    def create_dosage(self, token):
        
        n_ind = int(token.shape[0]/3)
        token = token.reshape(n_ind,3)

        if self.reference_allelle == 1:
            coeff = np.array([2,1,0])
            token_dosage = np.dot(token, coeff)

        elif self.reference_allelle == 2:
            coeff = np.array([0,1,2])
            token_dosage = np.dot(token, coeff)   
            self.dosage_relative_to = 4

        else:
            coeff = np.array([2,1,0])
            token_dosage = np.dot(token, coeff)
            
            minor_allele_freq = float(np.sum(token_dosage)/(2*n_ind))
            
            if minor_allele_freq < 0.5:
                return token_dosage
            else:
                coeff = np.array([0,1,2])
                token_dosage = np.dot(token, coeff)
                self.dosage_relative_to = 4 #change

        return token_dosage


    def process_file(self, filename, ch_nr, snp_stats_stem):
   
        filepath = os.path.join(self.input_directory, filename)

        utils.check_file_existence(filepath)

        ch_snp_stats_path = os.path.join(self.snp_stats_path, snp_stats_stem + str(ch_nr) + '.snp-stats.gz')
        utils.check_file_existence(ch_snp_stats_path)

        read = self.fileprocessingsystem.read_file(filename.split('.')[-1], filepath)
        length = len(read.readlines())
        read.seek(0)
        dosage_final = []


        for line in read:
            dosage_final.append(self.process_token(line, length))

        # Convert to numpy array because HDF accepts that
        self.map_file = np.asarray(self.map_file, dtype='S')
        np.save('./data/map_numpy.npy', self.map_file)
        dosage_final = np.asarray(dosage_final)
        np.save('./data/dosage_numpy.npy', dosage_final)
        print(dosage_final.shape)
        sys.stdout.write("\n")
        self.logger.warning('----DOSAGE FILE----')

        self.logger.warning('Dosage file for chromosome {} contains data for {} SNPs and {} patients'.format(ch_nr, dosage_final.shape[0], dosage_final.shape[1]))
        
        # Close the read
        read.close()
        
        snp_stats_data = self.read_snp_stats(ch_snp_stats_path)
        
        self.logger.warning('----MAP FILE----')
        if self.map_file.shape[0] != snp_stats_data.shape[0]:
            self.logger.warning('Map file shape: {}, snp_stats_data shape: {}'.format(self.map_file.shape, snp_stats_data.shape))
            self.logger.warning('!! Map file and snp_stats file do not have the same nr. of SNPs for Chromosome {}!'.format(ch_nr))
            sys.exit()
            
        else:
            self.map_file = np.column_stack((self.map_file,snp_stats_data[:,1:]))
            self.logger.warning('Map file shape: {}, snp_stats_data shape: {}'.format(self.map_file.shape, snp_stats_data.shape))
            self.logger.warning('Map file for chromosome {} contains {} variables for {} SNPs'.format(ch_nr, self.map_file.shape[1],self.map_file.shape[0]))

            
        return dosage_final
    
    
    def process_token(self, token, length):

        barLength = 100
        l = token.strip().split(' ')
        three_prob_data = l[self.header:]
        self.map_file.append(l[:self.header-2])
        
        data_to_convert = np.asarray(three_prob_data, dtype = float)

        count = 0
        if self.showProgress:
        
            count += 1
            progress = int(round(count/float(length)*barLength))
            if progress%1==0:
                sys.stdout.write("\r[{0}] {1} % ".format("#"*progress + "-"*(barLength-progress), progress))
                sys.stdout.flush()
        
        dosage = self.create_dosage(data_to_convert)
        #dosage_token.append(dosage)
        relative_al = l[self.dosage_relative_to]
        self.map_file[-1] += [relative_al,self.reference_allelle] #change

        return dosage.tolist()


    def read_snp_stats(self, filepath):
        '''
        Parts are hardcoded given that the column names are always in the same order
        '''
        snp_stats_data = []
        read = self.fileprocessingsystem.read_file(filepath.split('.')[-1], filepath)
 
        c = 0
        header = ['SNP_id','minor_allele', 'major_allele', 'MAF', 'HWE', 'information']
        next(read)
        for line in read:
            l = line.strip().split(' ')
            information = float(l[-1])
            hwe = float(l[-4])
            maf = float(l[-5])
            minor_al = str(l[6])
            major_al = str(l[7])
            snp_id = str(l[1])
            snp_stats_data.append([snp_id, minor_al, major_al, maf, hwe, information])
        snp_stats_data = np.asarray(snp_stats_data, dtype=object)

        read.close()
        return snp_stats_data
    
    def save_sample(self, path, filename, hdf_filepath):
        self.logger.warning('----SAMPLE FILE----')
        self.logger.warning('- Sample parametar = 1')
        if path[-1] == '/':
            filepath = path + filename
        else:
            filepath = path + '/' + filename
        
        read = self.fileprocessingsystem.read_file(filepath.split('.')[-1], filepath)

        header = read.readline().rstrip().split(' ')
        data = np.loadtxt(read, delimiter=' ', dtype=bytes).astype(str)
        #dt = np.dtype(dtype='|S125')
        data_sample = np.asarray(data[1:], dtype=object)
        read.close()
        dt = h5py.special_dtype(vlen=str)
       
        with h5py.File(hdf_filepath, 'a') as hdf:
            if '/Sample' in hdf.keys():
                self.logger.warning('!!Sample file is already stored in {} \n'.format(hdf_filepath))
                pass
            else:
                G = hdf.create_group('Sample')
            # G1 = hdf.create_group('Genetics')
                G.create_dataset('Sample', data=data_sample, dtype=dt, compression = 'lzf')
                G.attrs['sample_header'] = np.array(header, dtype=dt)
                self.logger.warning('Sample file saved. You can find it in {} under /Genetics/Sample \n'.format(hdf_filepath))
        #data_sample = data_sample.astype(np.str)
    
    
    def run_converter(self, filename, ch_nr, hdf_filepath, snp_stats_stem):
        hdf_group_name = 'Genetics/' + 'Ch' + ch_nr
        hdf_dataset_name = 'dosage' + ch_nr
        hdf_map_file_name = 'map' + ch_nr
        
        #first check if chromosome is already processed
        self.check_ch_processed(hdf_filepath, hdf_group_name)
        
        self.logger.warning('STATUS: Converting the data for Chromosome {}:'.format(ch_nr))
        start_time = time.time()
        dosage = self.process_file(filename, ch_nr, snp_stats_stem)
        self.logger.warning("STATUS: Data {} converted into dosages. Time needed:{} seconds ---\n".format(filename,time.time() - start_time)) 
        start_time = time.time()
   
        self.logger.warning('STATUS: Saving the data for Chromosome {}'.format(ch_nr))
        print(self.map_file)
        self.map_file = self.map_file.astype('S')
        print(self.map_file.dtype)
        self.save_file(hdf_filepath, hdf_group_name, hdf_dataset_name, hdf_map_file_name, dosage)
        self.logger.warning("STATUS: Data saved into file: {} under {}. Time needed: {} seconds ---".format(hdf_filepath, hdf_group_name, time.time()-start_time))


    def save_file(self, hdf_filepath, hdf_group_name, hdf_dataset_name, hdf_map_file_name, dosage):
        #dt = h5py.special_dtype(vlen=str)    
        #dt = h5py.vlen_dtype(np.dtype('S1'))
        with h5py.File(hdf_filepath, 'a') as hdf:
            G = hdf.create_group(hdf_group_name)
            G.create_dataset(hdf_map_file_name, data=self.map_file, dtype=h5py.string_dtype(), compression = 'lzf')
            #G.create_dataset(hdf_dataset_name, data=dosage, compression = 'lzf')


    def check_ch_processed(self, hdf_filepath, hdf_group_name):
        with h5py.File(hdf_filepath, 'a') as hdf:
            if hdf_group_name in hdf.keys(): #!!!!!!!ovdje staviti ako postoji u genetics
                self.logger.warning('!! There already exists group {} for the given Chromosome!'.format(hdf_group_name))
                sys.exit()
            else:
                pass



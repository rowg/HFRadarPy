import datetime as dt
import glob
import os
import re
import pandas as pd
import io
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from calc import dms2dd, createLonLatGridFromTopLeftPointWera, createLonLatGridFromBB
from pyproj import Geod
import math
import numpy as np

import logging

logger = logging.getLogger(__name__)


def divide_chunks(iterable, chunksize):
    """
    Yield successive n-sized chunks from a list

    Link: https://www.geeksforgeeks.org/break-list-chunks-size-n-python/

    Args:
        l (list): list to be broken down
        n (int): integer of the size chunks you want from the list

    Yields:
        list: a list containing lists of size, n
    """
    # looping till length l
    for i in range(0, len(iterable), chunksize):
        yield iterable[i: i + chunksize]


def list_files(types, main_dir, sub_directories=()):
    """
    Return a list of files given the directory of the files and extension type.
    You may also provide a list of sub_directories for the function to avoid.

    Args:
        types (str): file extension that you want to find
        main_dir (_type_): main directory that you want to recursively search for files
        sub_directories (tuple, optional):  Tuple containing strings of subdirectories you want to avoid. Defaults to ().

    Returns:
        list: list of files
    """
    file_list = []  # create empty list for finding files

    sub_dirs = [
        os.path.join(main_dir, o)
        for o in os.listdir(main_dir)
        if os.path.isdir(os.path.join(main_dir, o)) and o in sub_directories
    ]

    for sub in sub_dirs:
        for ext in types:
            file_list.extend(glob.glob(os.path.join(sub, ext)))
    file_list = sorted(file_list)
    return file_list


def list_to_dataframe(file_list):
    """
    Convert a list of ctf files, that are named in the the standard format 'year_month_day_hours' to a pandas dataframe

    Args:
        file_list (list): a list of files

    Returns:
        pd.DataFrame: Pandas DataFrame containing a list of files
    """
    df = pd.DataFrame(sorted(file_list), columns=["file"])
    try:
        df["time"] = df["file"].str.extract(r"(\d{4}_\d{2}_\d{2}_\d{4})")
        df["time"] = df["time"].apply(lambda x: dt.datetime.strptime(x, "%Y_%m_%d_%H%M"))
        df = df.set_index(["time"]).sort_index()
    except ValueError:
        logging.error("Cannot pass empty file_list to function. Returning empty dataframe.")
    return df


def timestamp_from_lluv_filename(filename):
    """
    Convert the string timestamp represented in CTF file names into a dt.datetime.

    Args:
        filename (str): filename

    Returns:
        dt.datetime: a datetime representation of the time included in the ctf filename
    """
    timestamp_regex = re.compile(r"\d{4}_\d{2}_\d{2}_\d{4}")
    mat_time = timestamp_regex.search(filename).group()
    timestamp = dt.datetime.strptime(mat_time, "%Y_%m_%d_%H%M")
    return timestamp

def addBoundingBoxMetadata(obj,lon_min,lon_max,lat_min,lat_max,grid_res=None):
    """
    This function adds metadata related to the bounding box to the input Radial or Total object.
    
    INPUTS:
        obj: Radial or Total object
        lon_min: minimum longitude of the bounding box
        lon_max: maximum longitude of the bounding box
        lat_min: minimum latitude of the bounding box
        lat_max: maximum latitude of the bounding box
        grid_res: grid resolution in km

        
    OUTPUTS:
        obj = Radial or Total object with metadata related to the bounding box
        
    """
    obj.metadata['BBminLongitude'] = str(lon_min) + ' deg'
    obj.metadata['BBmaxLongitude'] = str(lon_max) + ' deg'
    obj.metadata['BBminLatitude'] = str(lat_min) + ' deg'
    obj.metadata['BBmaxLatitude'] = str(lat_max) + ' deg'
    if grid_res:
        obj.metadata['GridSpacing'] = str(grid_res) + ' km'
    
    return obj
    
class fileParser(object):
    """
    A generic parser for the CODAR CTF, WERA CRAD and CUR file formats.
    """

    __metaclass__ = ABCMeta

    def __init__(self, fname=''):
        """
        Return an fileParser object
        """
        self.metadata = OrderedDict()
        self._tables = OrderedDict()
        if fname:
            split_path = os.path.split(fname)
            self.file_path = split_path[0]
            self.file_name = split_path[1]
            self.full_file = os.path.realpath(fname)
            extension = os.path.splitext(fname)[1]
            
            if (extension == '.ruv') or (extension == '.tuv'):
                self.CTFparser()
            elif extension == '.crad_ascii':
                self.CRADparser()
            elif extension == '.cur_asc':
                self.CURparser()  
    
    def CTFparser(self):
        """
        Return an fileParser object obtained by parsing CTF-LLUV files
        """
        table_count = 0
        table = False  # Set table to False. Once a table is found, switch to True.
        self.is_wera = False  # Default false. If 'WERA' is detected in the Manufacturer flag it is then set to True
        self.is_combined = False  # Default false. If radial combination is performed it is then set to True
        processing_info = []
        site_source = []

        with open(self.full_file, 'r', encoding='ISO-8859-1') as open_file:
            open_lluv = open_file.readlines()
            if any('%End:' in s or s.strip() == '%End' for s in open_lluv):  # if there is no %End: the file is corrupt!
                # Parse header and footer metadata
                for line in open_lluv:

                    # Fix for older WERA files
                    # Add a colon to the end of '%End'
                    if line.strip() == '%End':
                        line += ':'

                    if not table:  # If we are not looking at a table or a tables header information
                        if line.startswith('%%'):
                            if 'SiteSource' in line:
                                site_source.append(line)
                            else:
                                continue
                        elif line.startswith('%'):  # Parse the single commented header lines
                            key, value = self._parse_header_line(line)
                            if 'TableType' in line:  # Save this data as global header information
                                table = True  # we found a table
                                table_count = table_count + 1  # this is the nth table
                                table_data = u''
                                # self._data_header[table_count] = []
                                self._tables[str(table_count)] = OrderedDict()
                                self._tables[str(table_count)][key] = value
                                self._tables[str(table_count)]['_TableHeader'] = []
                            elif 'Manufacturer' in line:
                                if 'WERA' in value:
                                    self.is_wera = True
                                self.metadata[key] = value
                            elif 'SiteSource' in line:
                                site_source.append(value)
                            elif table_count > 0:
                                if key == 'ProcessingTool':
                                    processing_info.append(value)
                                else:
                                    self.metadata[key] = value
                            else:
                                self.metadata[key] = value
                    elif table:
                        if line.startswith(('%', ' %')):
                            if line.startswith(('%%', ' %%')):  # table header information
                                rep = {' comp': '_comp', ' Distance': '_Distance',' Ratio': '_Ratio',' (dB)': '_(dB)',' Width': '_Width', ' Resp': '_Resp', 'Value ': 'Value_','FOL ':'FOL_' }
                                rep = dict((re.escape(k), v) for k, v in rep.items())
                                pattern = re.compile('|'.join(rep.keys()))
                                temp = pattern.sub(lambda m: rep[re.escape(m.group(0))], line).strip('% \n')
                                temp = [x.replace('_', ' ') for x in re.sub(' +', ' ', temp).split(' ')]  # Get rid of underscores
                                # temp[0] = '%%   {}'.format(temp[0])

                                self._tables[str(table_count)]['_TableHeader'].append(temp)
                            else:  # Table metadata and diagnostic data are prepended by at least 1 % sign
                                if len(line.split(':')) == 1:  # Diagnostic Data
                                    line = line.replace('%', '').strip()
                                    table_data += '{}\n'.format(line)
                                else:  # Table data
                                    key, value = self._parse_header_line(line)
                                    # if 'TableColumnTypes' not in self._tables[str(table_count)]:
                                    #     raise ValueError("TableColumnTypes not defined")
                                    if 'TableEnd' in line:
                                        if 'TableColumnTypes' in self._tables[str(table_count)]:
                                            # use pandas read_csv because it interprets the datatype for each column of the csv
                                            tdf = pd.read_csv(
                                                io.StringIO(table_data),
                                                sep=' ',
                                                header=None,
                                                names=self._tables[str(table_count)]['TableColumnTypes'].split(),
                                                skipinitialspace=True
                                            )
                                        else:
                                            tdf = pd.DataFrame()

                                        self._tables[str(table_count)]['data'] = tdf
                                        table = False
                                    else:
                                        key, value = self._parse_header_line(line)
                                        self._tables[str(table_count)][key] = value
                        else:  # Uncommented lines are the main data table.
                            table_data += '{}'.format(line)
                self.metadata['ProcessingTool'] = processing_info
                if site_source:
                    self.metadata['SiteSource'] = site_source
                # if self.is_wera:
                #     self._tables['1']['data'] = pd.read_csv(io.StringIO(table_data),
                #                                             sep=' ',
                #                                             header=None,
                #                                             names=['LOND', 'LATD', 'VELU', 'VELV', 'VFLG', 'EACC', 'RNGE', 'BEAR', 'VELO', 'HEAD'],  # WERA has incorrect TableColumnTypes in their files.....
                #                                             skipinitialspace=True, )
                self._iscorrupt = False
            else:
                logging.error('{}: File corrupt. Skipping to next file.'.format(self.full_file))
                self._iscorrupt = True
        try:
            self.time = dt.datetime(*[int(s) for s in self.metadata['TimeStamp'].split()])
        except KeyError:
            pass
        
    def CRADparser(self):
        """
        Return an fileParser object obtained by parsing WERA CRAD radial files
        """
        # Load the WERA crad_ascii Data with this generic CRAD parsing routine below
        table_count = 0
        table = False  # Set table to False. Once a table is found, switch to True.
        self.is_wera = True  # Default True
        processing_info = []
        site_source = []

        with open(self.full_file, 'r') as open_file:
            open_crad = open_file.readlines()
            open_crad = [i.lstrip() for i in open_crad]
            # Parse header
            header= str(''.join(open_crad[0 : 9])).replace("\n", " ").strip()
            self.metadata = self._parse_crad_header(header)
            # Read data content
            table = True  # we found a table
            table_count = table_count + 1  # this is the nth table
            table_data = u''
            self._tables[str(table_count)] = OrderedDict()
            self._tables[str(table_count)]['TableType'] = 'CRAD'
            self._tables[str(table_count)]['_TableHeader'] = ['Top-Left gridpoint Latitude','Top-Left gridpoint Longitude','Number of measurements','Average (SignalToNoise * RadialVelocity)','Average (SignalToNoise * SquaredRadialVelocity)','Sum over all Signal-to-Nois-Ratio','Overall Power of the gridcell']
            self._tables[str(table_count)]['TableColumnTypes'] = 'LatC LonC KUR SNV SNS SNR PWR'
            table_data = ''.join(open_crad[9:])
            tdf = pd.read_csv(
                io.StringIO(table_data),
                sep=' ',
                header=None,
                names=self._tables[str(table_count)]['TableColumnTypes'].split(),
                skipinitialspace=True
            )
            self._tables[str(table_count)]['data'] = tdf
            
        # Get the indexes of rows for which Kur is 0 (i.e. no measurements)
        indexNames = tdf[ tdf['KUR'] == 0 ].index
        # Delete these row indexes from DataFrame
        tdf.drop(indexNames , inplace=True)
        tdf.reset_index(level=None, drop=False, inplace=True)
        
        # Get the site coordinates
        siteLon = dms2dd(list(map(int,self.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
        siteLat = dms2dd(list(map(int,self.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][:-2].split('-'))))
        if self.metadata['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'S':
            siteLat = -siteLat
        if self.metadata['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'][-1] == 'W':
            siteLon = -siteLon
        
        # Use WGS84 ellipsoid
        g = Geod(ellps='WGS84')
        self.metadata['GreatCircle'] = '"WGS84"' + ' ' + str(g.a) + '  ' + str(1/g.f)
        
        # Parse data content
        radialData = tdf.apply(lambda x: self._parse_crad_data(x,siteLon,siteLat,g), axis=1)
        # Assign column names to the combination DataFrame
        radialData.columns = ['LOND','LATD','VELU','VELV','VELO','HEAD','HCSS','EACC']   
        table = True  # we found a table
        table_count = table_count + 1  # this is the nth table
        table_data = u''
        self._tables[str(table_count)] = OrderedDict()
        self._tables[str(table_count)]['TableType'] = 'LLUV'
        self._tables[str(table_count)]['_TableHeader'] = ['Longitude', 'Latitude', 'U comp', 'V comp', 'Velocity', 'Direction', 'Variance', 'Accuracy']
        self._tables[str(table_count)]['TableColumnTypes'] = 'LOND','LATD','VELU','VELV','VELO','HEAD','HCSS','EACC'
        self._tables[str(table_count)]['data'] = radialData
        
        self._iscorrupt = False
    
        try:
            self.time = dt.datetime.strptime(self.metadata['DateOfMeasurement'], '%d-%b-%y %H:%M %Z')
        except KeyError:
            pass
        
    def CURparser(self):
        """
        Return an fileParser object obtained by parsing WERA CUR total files
        """
        table_count = 0
        table = False  # Set table to False. Once a table is found, switch to True.
        self.is_wera = True  # Default True
        self.is_combined = False  # Default false. If radial combination is performed it is then set to True
        processing_info = []
        site_source = []

        with open(self.full_file, 'r') as open_file:
            open_crad = open_file.readlines()
            open_crad = [i.lstrip() for i in open_crad]
            # Parse header
            numStation = int(open_crad[0])
            header= open_crad[0 : numStation+7]
            header = [elem for elem in header if elem.strip()]
            self.metadata,self.site_source = self._parse_cur_header(header)
            # Read data content
            table = True  # we found a table
            table_count = table_count + 1  # this is the nth table
            table_data = u''
            self._tables[str(table_count)] = OrderedDict()
            self._tables[str(table_count)]['TableType'] = 'CUR'
            self._tables[str(table_count)]['_TableHeader'] = ['X-direction Index','Y-direction Index','Velocity U component [m/s]','Velocity V component [m/s]','Klass','Accuracy of U component [m/s]','Accuracy of V component [m/s]']
            self._tables[str(table_count)]['TableColumnTypes'] = 'IX IY U V KL Acc_U Acc_V'
            table_data = ''.join(open_crad[numStation+9:])
            tdf = pd.read_csv(
                io.StringIO(table_data),
                sep=' ',
                header=None,
                names=self._tables[str(table_count)]['TableColumnTypes'].split(),
                skipinitialspace=True
            )
            self._tables[str(table_count)]['data'] = tdf
            
        # Get grid information
        cellSize = float(self.metadata['DGT'].split()[0])
        topLeftLon = float(self.metadata['TopLeftLongitude'])
        topLeftLat = float(self.metadata['TopLeftLatitude'])
        nx = int(self.metadata['NX'])
        ny = int(self.metadata['NY'])
        
        # Generate grid coordinates
        gridGS = createLonLatGridFromTopLeftPointWera(topLeftLon, topLeftLat, cellSize, nx, ny)
        # extract longitudes and latitude from grid GeoSeries and insert them into numpy arrays
        lonVec = np.unique(gridGS.x.to_numpy())
        latVec = np.flipud(np.unique(gridGS.y.to_numpy()))
        # manage antimeridian crossing
        lonVec = np.concatenate((lonVec[lonVec>=0],lonVec[lonVec<0]))
                
        # Use WGS84 ellipsoid
        g = Geod(ellps='WGS84')
        self.metadata['GreatCircle'] = '"WGS84"' + ' ' + str(g.a) + '  ' + str(1/g.f)
            
        # Parse data content
        totalData = tdf.apply(lambda x: self._parse_cur_data(x,lonVec,latVec), axis=1)
        # Assign column names to the combination DataFrame
        totalData.columns = ['LOND','LATD','VELU','VELV','VELO','HEAD','UACC','VACC']   
        table = True  # we found a table
        table_count = table_count + 1  # this is the nth table
        table_data = u''
        self._tables[str(table_count)] = OrderedDict()
        self._tables[str(table_count)]['TableType'] = 'LLUV'
        self._tables[str(table_count)]['_TableHeader'] = ['Longitude', 'Latitude', 'U comp', 'V comp', 'Velocity', 'Heading', 'U accuracy', 'V accuracy']
        self._tables[str(table_count)]['TableColumnTypes'] = 'LOND','LATD','VELU','VELV','VELO','HEAD','UACC','VACC'
        self._tables[str(table_count)]['data'] = totalData
        
        self._iscorrupt = False
    
        try:
            self.time = dt.datetime.strptime(self.site_source['DateOfMeasurement'][0], '%d-%b-%Y %H:%M %Z')
        except KeyError:
            pass

    def is_valid(self, table='1'):
        """
        Check if the data table for the file contains data

        Args:
            table (str, optional): string containing the table number to validate. Defaults to '1'.

        Returns:
            bool: True or False if data is present
        """
        try:
            return not self._tables[table]['data'].empty
        except:
            return False

    @staticmethod
    def _parse_header_line(line):
        """
        Parse a line into a key, value
        
        INPUT:
            line: line from a text file
        
        OUTPUT:
            key,value: tuple containing the key, value for the line
        """

        line = line.replace('%', '')  # Strip the % sign from the line
        line = line.replace('\n', '')  # Strip the new line character from the end of the line
        line_split = line.split(':')
        key = line_split[0]  # save key variable
        value = line_split[1].strip()  # save value variable and strip whitespace from value
        return key, value
    
    @staticmethod
    def _parse_crad_header(header):
        """
        Parse CRAD header into key, value pairs and assigns them to the metadata field
        
        INPUT:
            header: string containing the crad header
        
        OUTPUT:
            metadataDict: dictionary containing the key, value pairs from the header
        """

        # create output Dictionary
        metadataDict =  OrderedDict()
        
        # Split header into word list (blank space, comma, column and equal as separators)
        header = header.replace('RA TE', 'RATE')
        wordList=re.split(r"\s+|,+|:+|=",header)
        wordList = list(filter(None, wordList))     # Remove empty elements

        # Parse header
        if 'SAMPLES' in wordList:
            metadataDict['Samples'] = wordList[wordList.index("SAMPLES")-1].strip()
            metadataDict['DateOfMeasurement'] = wordList[wordList.index("SAMPLES")+1].strip() + ' ' + \
                                                wordList[wordList.index("SAMPLES")+2].strip() + ':' + \
                                                wordList[wordList.index("SAMPLES")+3].strip() + ' ' + \
                                                wordList[wordList.index("SAMPLES")+4].strip()
            metadataDict['TimeZone'] = metadataDict['DateOfMeasurement'].split()[2]
            metadataDict['StationName'] = wordList[wordList.index("SAMPLES")+5].strip()
        if 'FREQUENZ' in wordList:
            metadataDict['FileType'] = wordList[wordList.index("FREQUENZ")-1].strip()
            metadataDict['CenterFrequency'] = wordList[wordList.index("FREQUENZ")+1].strip()
        if 'YEAR' in wordList:
            metadataDict['Year'] = wordList[wordList.index("YEAR")+1].strip()
        if 'RANGE' in wordList:
            metadataDict['Range'] = wordList[wordList.index("RANGE")+1].strip() + ' ' + \
                                    wordList[wordList.index("RANGE")+2].strip()
        if 'TRUENORTH' in wordList:
            metadataDict['TrueNorth'] = wordList[wordList.index("TRUENORTH")+1].strip() + ' ' + \
                                        wordList[wordList.index("TRUENORTH")+2].strip()
        if 'RATE' in wordList:
            metadataDict['ChirpRate'] = wordList[wordList.index("RATE")+1].strip()
        if 'NRRANGES' in wordList:
            metadataDict['NumberOfRangeCells'] = wordList[wordList.index("NRRANGES")+1].strip()
            metadataDict['NumberOfCoherentUSORTfilesSinceMeasurementStart'] = wordList[wordList.index("NRRANGES")+2].strip()
        if any("BREITE" in item for item in wordList):
            idx = wordList.index([item for item in wordList if "BREITE" in item][0])
            metadataDict['Latitude(deg-min-sec)OfTheCenterOfTheReceiveArray'] = wordList[idx+1].strip() + \
                                                                                wordList[idx+2].strip() + ' ' + \
                                                                                wordList[idx+3].strip()
            if 'LAENGE' in wordList:  
                if 'EBREITE' in wordList:
                    metadataDict['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'] = wordList[wordList.index("LAENGE")+1].strip() + ' E'
                elif 'WBREITE' in wordList:
                    metadataDict['Longitude(deg-min-sec)OfTheCenterOfTheReceiveArray'] = wordList[wordList.index("LAENGE")+1].strip() + ' W'
        if 'NTL' in wordList:
            metadataDict['NTL'] = wordList[wordList.index("NTL")+1].strip()
        if 'NFTL' in wordList:
            metadataDict['NFTL'] = wordList[wordList.index("NFTL")+1].strip()
        if 'nx' in wordList:
            metadataDict['nx'] = wordList[wordList.index("nx")+1].strip()
        if 'ny' in wordList:
            metadataDict['ny'] = wordList[wordList.index("ny")+1].strip()
        if 'OFFSET' in wordList:
            metadataDict['Offset'] = wordList[wordList.index("OFFSET")+1].replace('RXOFFSET','')
        if any("RXOFFSET" in item for item in wordList):
            idx = wordList.index([item for item in wordList if "RXOFFSET" in item][0])
            metadataDict['RXOffset'] = wordList[idx+1].replace('SS','')
        if any("SS" in item for item in wordList):
            idx = wordList.index([item for item in wordList if "SS" in item][0])
            metadataDict['SS'] = wordList[idx+1].strip()
        if 'HD' in wordList:
            idx1 = wordList.index("HD")+1
            if 'RFI2N' in wordList:
                idx2 = wordList.index("RFI2N")
            elif 'NCOV' in wordList:
                idx2 = wordList.index("NCOV")
            metadataDict['HD'] = ''
            for idx in range(idx1,idx2):
                metadataDict['HD'] = metadataDict['HD'] + wordList[idx] + ' '
            metadataDict['HD'].strip()
        if 'RFI2N' in wordList:
            metadataDict['RFI2N'] = wordList[wordList.index("RFI2N")+1].strip()
        if 'NCOV' in wordList:
            metadataDict['NCOV'] = wordList[wordList.index("NCOV")+1].strip()
        if 'LAT' in wordList:
            metadataDict['TopLeftLatitude'] = wordList[wordList.index("LAT")+1].strip()
        if 'LON' in wordList:
            metadataDict['TopLeftLongitude'] = wordList[wordList.index("LON")+1].strip()
        if 'DGT' in wordList:
            metadataDict['DGT'] = wordList[wordList.index("DGT")+1].strip()
        metadataDict['NumberOfSeries'] = wordList[-6]
        metadataDict['NumberOfAntennas'] = wordList[-5].replace('-','')
        metadataDict['CartesianRadials'] = wordList[-4].replace('-','')      
        
        return metadataDict
    
    @staticmethod
    def _parse_crad_data(cellData,siteLon,siteLat,g):
        """
        Parse crad data content into Radial parameters (i.e. CTF-like)
        
        INPUT:
            cellData: aeries containing the cell data
            siteLon: longitude of the radar site
            siteLat: latitude of the radar site
            g: Geod object with CRS
            
        OUTPUT:
            radialData: Series containing the Radial parameters (i.e. CTF-like)
        """

        # create output Series
        radialData = pd.Series(np.nan,index=range(8))
        
        # Parse data
        radialData.loc[0] = np.rad2deg(cellData['LonC'])
        radialData.loc[1] = np.rad2deg(cellData['LatC'])
        radialData.loc[4] = cellData['SNV'] / cellData['SNR']
        radialData.loc[5],az21,dist = g.inv(siteLon,siteLat,radialData.loc[0],radialData.loc[1])
        if radialData.loc[5] <0:
            radialData.loc[5] += 360    # keep angles clockwise from true North
        radialData.loc[2] = radialData.loc[4] * math.sin(math.radians(radialData.loc[5]))
        radialData.loc[3] = radialData.loc[4] * math.cos(math.radians(radialData.loc[5]))
        radialData.loc[6] = cellData['SNS'] / cellData['SNR']
        radialData.loc[7] = radialData.loc[6] / math.sqrt(cellData['KUR'])              
        
        return radialData
    
    @staticmethod
    def _parse_cur_header(header):
        """
        Parse CUR header into key, value pairs and assigns them to the metadata field
        
        INPUT:
            header: list containing the cur header
        
        OUTPUT:
            metadataDict: dictionary containing the key, value pairs from the header
            site_source: DataFrame containing the site source information
        """

        # create output Dictionary
        metadataDict =  OrderedDict()

        # Parse first part of the header (contributing radial sites)
        metadataDict['NumberOfContributingRadials'] = header[0]
        lineNumber = int(metadataDict['NumberOfContributingRadials'])
        metadataDict['SiteSource'] = ['DateOfMeasurement Name Lat Lon Coverage(s)']
        for i in range(lineNumber):
            metadataDict['SiteSource'].append(header[i+1])
        
        # Parse second part of the header (regular structure)
        lineNumber += 3
        line = header[lineNumber].strip()
        elements = line.split()
        metadataDict['TopLeftLatitude'] = elements[0]
        metadataDict['TopLeftLongitude'] = elements[1]
        metadataDict['DGT'] = elements[2]
        metadataDict['NX'] = elements[3]
        metadataDict['NY'] = elements[4]
        lineNumber += 1
        metadataDict['NumberOfEntries'] = header[lineNumber]   
        
        # create output DataFrame
        site_source = pd.DataFrame(np.nan,index=range(len(metadataDict['SiteSource'])-1),columns=['#', 'Name', 'Lat', 'Lon', 'Coverage(s)', 'DateOfMeasurement', 'RngStep(km)', 'Pattern', 'AntBearing(NCW)'])
        
        # Parse site source information
        dateOfMeasurement = []
        name = []
        siteLat = []
        siteLon = []
        coverage = []
        
        for ss in metadataDict['SiteSource']:
            if not 'DateOfMeasurement' in ss:
                dateOfMeasurement.append(ss.split()[0] + ' ' + ss.split()[1] + ' ' + ss.split()[2])
                name.append(ss.split()[3])
                if ss.split()[5] == 'North':
                    siteLat.append(float(ss.split()[4]))
                elif ss.split()[5] == 'South':
                    siteLon.append(-float(ss.split()[4]))
                if ss.split()[7] == 'East':
                    siteLon.append(float(ss.split()[6]))
                elif ss.split()[7] == 'West':
                    siteLon.append(-float(ss.split()[6]))
                coverage.append(ss.split()[10])
                
        site_source['#'] = np.arange(len(site_source.index))+1
        site_source['Name'] = name
        site_source['Lat'] = siteLat
        site_source['Lon'] = siteLon
        site_source['Coverage(s)'] = coverage
        site_source['DateOfMeasurement'] = dateOfMeasurement
        
        return metadataDict, site_source
    
    @staticmethod
    def _parse_cur_data(cellData,lonVec,latVec):
        """
        Parse cur data content into Total parameters (i.e. CTF-like)
        
        INPUT:
            cellData: Series containing the cell data
            lonVec: array containing the longitudes of the geographical grid
            latVec: array containing the latitudes of the geographical grid
            
        OUTPUT:
            totalData: Series containing the Total parameters (i.e. CTF-like)
        """

        # create output Series
        totalData = pd.Series(np.nan,index=range(7))
        
        # Parse data
        totalData.loc[0] = lonVec[int(cellData['IX']-1)]
        totalData.loc[1] = latVec[int(cellData['IY']-1)]
        totalData.loc[2] = cellData['U']
        totalData.loc[3] = cellData['V']
        totalData.loc[4] = np.sqrt(cellData['U']**2 + cellData['V']**2)
        totalData.loc[5] = (360 + np.arctan2(cellData['U'],cellData['V']) *180 / np.pi) % 360
        totalData.loc[6] = cellData['Acc_U']
        totalData.loc[7] = cellData['Acc_V']
        
        return totalData

    @abstractmethod
    def file_type(self):
        """Return a string representing the type of file this is."""
        return self.metadata['FileType']

    def replace_invalid_values(self, values=[999.00, 1080.0]):
        """
        Convert invalid CODAR values to NaN
        
        INPUT:
            df: dataframe
            values: list of CODAR fill values that reflect non calculable values
            
        OUTPUT:
            dataframe with invalid values set to NaN
        """
        logging.info('Replacing invalid values {} with NaN'.format(values))
        self.data.replace(values, np.nan, inplace=True)

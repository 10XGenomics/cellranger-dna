__doc__="""Set of methods that override or augment tenkit.bio_io"""
from crdna.constants import PROCESSED_BARCODE_TAG, RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG

def get_read_raw_barcode(read):
    try:
        return read.opt(RAW_BARCODE_TAG)
    except KeyError:
        return None

def get_read_barcode_qual(read):
    try:
        return read.opt(RAW_BARCODE_QUAL_TAG)
    except KeyError:
        return None

def read_has_barcode(read):
    try:
        read.opt(RAW_BARCODE_TAG)
        return True
    except KeyError:
        return False

def get_read_barcode(read):
    '''Get the 10X barcode sequence for a read.  Returns None if no barcode attached
       or a non-whitelist sequence was observed '''
    try:
        r = read.opt(PROCESSED_BARCODE_TAG)
        if r == '':
            return None
        else:
            return r
    except KeyError:
        return None

def get_read_barcode_or_raw(read):
    try:
        r = read.opt(PROCESSED_BARCODE_TAG)
        if r == '':
            r = read.opt(RAW_BARCODE_TAG)
            return r
        else:
            return r
    except:
        r = read.opt(RAW_BARCODE_TAG)
        return r

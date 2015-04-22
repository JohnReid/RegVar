class Config(object):
    DEBUG = False
    TESTING = False
    UCSC_DIR = '/home/john/Data/UCSC'
    MAFFILES = (
        '/home/john/Data/UCSC/goldenPath/hg38/multiz20way/maf/chrY.maf.bz2',
    )
    CHOP = False


class ProductionConfig(Config):
    DATABASE_URI = 'mysql://user@localhost/foo'


class DevelopmentConfig(Config):
    DEBUG = True


class TestingConfig(Config):
    TESTING = True

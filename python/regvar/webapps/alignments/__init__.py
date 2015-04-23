class Config(object):
    DEBUG = False
    TESTING = False
    UCSC_DIR = '/home/john/Data/UCSC'
    HOST = '0.0.0.0'
    PORT = 9083


class ProductionConfig(Config):
    pass


class DevelopmentConfig(Config):
    DEBUG = True


class TestingConfig(Config):
    TESTING = True

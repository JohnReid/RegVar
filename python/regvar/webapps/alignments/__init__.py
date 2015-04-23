class Config(object):
    DEBUG = False
    TESTING = False
    UCSC_DIR = '/home/john/Data/UCSC'


class ProductionConfig(Config):
    pass


class DevelopmentConfig(Config):
    DEBUG = True


class TestingConfig(Config):
    TESTING = True

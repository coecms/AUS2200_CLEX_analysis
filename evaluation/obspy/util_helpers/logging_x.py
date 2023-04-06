import logging as pylog

def_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

logger_levels = { 'CRITICAL': pylog.CRITICAL,  # 50
                  'DEBUG': pylog.DEBUG,  # 10
                  'ERROR': pylog.ERROR,  # 40
                  'FATAL': pylog.FATAL,  # 50
                  'INFO': pylog.INFO,    # 20
                  'WARN': pylog.WARN,    # 30
                  'WARNING': pylog.WARNING  } 

def getLogger(id, level=pylog.WARN, fmt=def_fmt):
    """
    Make it easy to set up a logger
    """

    logger = pylog.getLogger(id)
    if isinstance(level,str):
        level = logger_levels[ level.upper() ]

    logger.setLevel(level)

    sh = pylog.StreamHandler()
    sh.setLevel(pylog.DEBUG)    # defined for all levels
    formatter = pylog.Formatter(fmt)
    sh.setFormatter(formatter)

    logger.addHandler(sh)

    return logger

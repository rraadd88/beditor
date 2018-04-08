import logging

def get_logger(argv=None):
    """
    Initiates logging information in a pre-defined format. 
    """
    import logging
    import datetime
    log_format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s'
    #'[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s'

    logging.basicConfig(format=log_format,
                        level=logging.DEBUG,)
    if not argv is None:
        log_fh="%s_%s" % (make_pathable_string(str(datetime.datetime.now())),'_'.join(argv).replace('/','|'))
        print(log_fh)
        logging.basicConfig(filename=log_fh)
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter(log_format)
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
    # logging.info('#START')
    return logging
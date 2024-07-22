from . import partitionModel_RECON2, partitionModel_Human1

import logging
logger = logging.getLogger('compass')

def partition_model(args):

    if args['model'] == 'RECON2_mat':
        return partitionModel_RECON2.partition_model(args)
    elif args['model'] == 'Human1':
        return partitionModel_Human1.partition_model(args)
    else:
        model = args['model']
        logger.info(f'{model} is not supported by Module-Compass')
        return (None, None)
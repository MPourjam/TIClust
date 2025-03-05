import logging


ticlust_logger = logging.getLogger(__name__)
ticlust_logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    '%(asctime)s - %(name)s:'
    '%(filename)s:%(lineno)d'
    ' - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
handler.setFormatter(formatter)
ticlust_logger.addHandler(handler)

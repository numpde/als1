# RA, 2021-03-17

from bugs import *
from twig import log
import tensorflow as tf

# ragged matmul
# https://github.com/tensorflow/tensorflow/issues/28109

# Note that `ragged_rank=1`

rt = tf.ragged.constant(
    [
        [
            [1, 2, 3],
            [4, 5, 6]
        ],
        [[1, 2, 3]]
    ],
    ragged_rank=1
)

t = tf.constant([[1], [2], [3]])

log.info(rt)
log.info(t)

r = tf.ragged.map_flat_values(tf.matmul, rt, t)
log.info(r)

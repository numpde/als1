# RA, 2021-03-19

import pandas as pd
import tensorflow as tf
from twig import log

# tf.get_logger().setLevel('ERROR')  # doesn't help

ds1 = tf.data.Dataset.from_tensor_slices([1, 2, 3, 4])
# .shuffle(buffer_size=10, seed=43, reshuffle_each_iteration=False)
ds2 = tf.data.Dataset.from_tensor_slices(list("ABCD"))
# .shuffle(buffer_size=10, seed=43, reshuffle_each_iteration=False)

log.info(list(ds1.take(4).as_numpy_iterator()))
log.info(list(ds2.take(4).as_numpy_iterator()))

log.info(list((tf.data.Dataset.zip((ds1, ds1))).take(4).as_numpy_iterator()))
log.info(list((tf.data.Dataset.zip((ds1, ds2))).take(4).as_numpy_iterator()))

s = pd.Series([6, 7, 8, 9])
# can use `lambda: iter(s)` or `s.__iter__`
ds3 = tf.data.Dataset.from_generator(s.__iter__, output_signature=tf.TensorSpec(shape=None, dtype=tf.int32))
# Note: no overflow
log.info(list(ds3.take(6).as_numpy_iterator()))
# With repeat:
log.info(list(ds3.repeat(2).take(6).as_numpy_iterator()))

df = pd.DataFrame(data={'x': [10, 20, 30], 'y': [40, 50, 60]})
ds4 = tf.data.Dataset.from_generator(df.values.__iter__, output_signature=tf.TensorSpec(shape=(2,), dtype=tf.float32))
log.info(list(ds4.take(6).as_numpy_iterator()))

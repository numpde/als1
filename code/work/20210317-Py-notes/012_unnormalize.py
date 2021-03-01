# RA, 2021-04-11

import tensorflow as tf

x = tf.constant([[1, 2, 3], [2, 4, 6], [1, 4, 8]], dtype=tf.float32)

n = tf.reduce_sum(x, axis=1, keepdims=True)

y = x / n

print(all(tf.reduce_sum(y, axis=1) == 1))

y = y * n

print(float(tf.reduce_max(tf.abs(x - y))))

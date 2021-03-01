# RA, 2021-04-11

# ?

import tensorflow as tf

a = tf.constant([[1] * 4, [2] * 4])
b = tf.constant([[3], [4]])

tf.RaggedTensor.from_tensor([a, b])

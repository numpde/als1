# RA, 2021-03-17

import tensorflow as tf

t = tf.constant([[1, 2, 3], [4, 5, 6]])

tt = tf.split(t, axis=1, num_or_size_splits=t.shape[1])

print("Original shape:", t.shape)
print("Split shapes:", [t.shape for t in tt])

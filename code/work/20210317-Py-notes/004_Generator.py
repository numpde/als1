# RA, 2021-03-19


import tensorflow as tf

data = list("ABCDEF")


class G:
    signature = (tf.TensorSpec(shape=(), dtype=tf.int32), tf.TensorSpec(shape=(1, None), dtype=tf.int32))

    def __init__(self):
        self.ncalls = 0

    def __call__(self):
        while True:
            self.ncalls += 1
            yield (self.ncalls, [[self.ncalls] * self.ncalls])


ds = tf.data.Dataset.from_generator(G(), output_signature=G.signature)

print(
    *list(ds.take(3)),
    sep='\n'
)

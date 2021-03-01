# RA, 2021-03-17

"""
Making a custom keras layer that does not complain.

The following Variables were used a Lambda layer's call (tf.linalg.matmul_549), but
are not present in its tracked objects:
  <tf.Variable 'fex_GO:2000979:0' shape=(1, 1) dtype=float32>
It is possible that this is intended behavior, but it is more likely
an omission. This is a strong indication that this layer should be
formulated as a subclassed Layer rather than a Lambda layer.

Solution: use `call` instead of `__call__`.
"""

from bugs import *
from twig import log
import tensorflow as tf

n_samples = 6
n_genes = 3

df = np.random.random_integers(low=1, high=9, size=[n_samples, n_genes])
df = tf.convert_to_tensor(df, dtype='float32')
ds = tf.data.Dataset.from_tensor_slices(df).batch(batch_size=2).prefetch(1)

class Layer(tf.keras.layers.Layer):
    def __init__(self):
        super().__init__()

    def build(self, input_shape):
        log.info(f"build with shape = {input_shape}")
        self.m = self.add_weight("m", shape=(input_shape[-1], input_shape[-1]), trainable=True, initializer='uniform')

    def call(self, x, training=False):
        return tf.matmul(x, self.m)


layer = Layer()

# This does not generate the warning
(layer(first(ds)))

model = tf.keras.Sequential([Layer()])

opt = tf.keras.optimizers.Adagrad(learning_rate=1)
model.compile(optimizer=opt, loss=tf.keras.losses.MeanSquaredError())

# This should learn the Identity matrix
model.fit(x=ds.repeat(100).map(lambda x: (x, x)), epochs=2)

print(layer.get_weights())
print(model.get_weights())

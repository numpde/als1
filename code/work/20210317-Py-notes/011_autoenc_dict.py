# RA, 2021-04-11

"""
Autoencoder with `dict` intermediate.
"""

import pandas as pd
import numpy as np
import tensorflow as tf
import contextlib


class Encoder(tf.keras.layers.Layer):
    def __init__(self):
        super(Encoder, self).__init__()
        self.w = self.add_weight(name="w", initializer=tf.keras.initializers.Zeros)

    def call(self, inputs, training=None, mask=None):
        inputs = inputs * self.w
        return dict(value=inputs, info="hello worlds.")


class Decoder(tf.keras.layers.Layer):
    def call(self, inputs, training=None, mask=None):
        assert isinstance(inputs, dict)
        return inputs['value']


input = tf.keras.Input(shape=())
encoder = Encoder()
decoder = Decoder()
model = tf.keras.Model(inputs=input, outputs=decoder(encoder(input)))

model(911)

with contextlib.suppress(AttributeError):
    # This doesn't work
    model.call(911)

w = 2

X = pd.DataFrame({'x': [1, 2, 3, 4]})
y = X.x * 2

model.compile(loss="mse", optimizer="sgd")
model.fit(X, y, epochs=10, verbose=0, batch_size=1)

assert np.isclose(w, encoder.w.numpy())


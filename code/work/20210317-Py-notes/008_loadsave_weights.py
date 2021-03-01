# RA, 2021-03-28

"""
model.save_weights
model.load_weights
mystery
"""

from tempfile import NamedTemporaryFile
import tensorflow as tf
import numpy as np
import h5py


class M(tf.keras.Model):
    def __init__(self):
        super().__init__()
        self.w = self.add_weight(name='w', shape=())

    def call(self, x, *arg):
        return self.w.numpy()


m = M()

with NamedTemporaryFile(mode='w') as f:
    m.save_weights(f.name, save_format='h5')
    print(m(1), m.get_weights())
    # 0.14819694 [0.14819694]

    m.set_weights(np.array([43]))
    print(m(1), m.get_weights())
    # 43.0 [43.0]

    m.load_weights(f.name)
    print(m(1), m.get_weights())
    # 43.0 [43.0]

    with h5py.File(f.name, mode='r') as h5:
        assert 0 == len(list(h5))

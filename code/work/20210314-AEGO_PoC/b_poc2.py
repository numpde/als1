# RA, 2021-03-16

# # # WIP # # #

import inspect

from bugs import *
from twig import log
from plox import rcParam
from progressbar import progressbar
from tcga.utils import First, Peek, Now
from networkx.algorithms.traversal.depth_first_search import dfs_tree

import obonet
import networkx as nx

if __name__ == '__main__':
    # Path for checkpoints etc
    logpath = mkdir(Path(__file__).with_suffix('') / f"run/{Now()}")

    # Copy of myself
    with (logpath / f"copy__{Path(__file__).name}").open(mode='w') as fd:
        fd.write(inspect.getsource(inspect.getmodule(inspect.stack()[0].frame)))


class Setup:
    URLS = {
        'expr': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/matrix.csv",
        'meta': "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hip_10x/metadata.csv",
        'obo': "http://purl.obolibrary.org/obo/go.obo",
    }

    # http://www.informatics.jax.org/vocab/gene_ontology/GO:0007417
    # central nervous system development (~900 genes)
    root = 'GO:0007417'

    # Fraction of samples to keep
    frac_samples = 0.03

    glob = Path(__file__).parent.glob

    sym2cat = pd.read_csv(unlist1(glob("../../data/*/*/goa_mouse_sym2cat.txt.gz")), sep='\t')

    from tcga.utils import download
    download = download.to(abs_path=unlist1(glob("../../data/*/Mouse-WCH*/*/download_cache")))

    with download(URLS['obo']).now.open() as fd:
        G = obonet.read_obo(fd)
        assert isinstance(G, nx.MultiDiGraph)

    # Take the subtree rooted at the `root` term
    G = G.subgraph(dfs_tree(G.reverse(), root))

    # Keep only the genes implicated in those GO terms
    sym2cat = sym2cat[sym2cat.goid.isin(G.nodes)]
    genes = sorted(set(sym2cat.symbol))

    path_to_reduced = mkdir(Path(__file__).with_suffix('') / "reduced_dataset")

    @classmethod
    def get_incidence_matrix(cls, genes):
        # Genes x GOID incidence matrix
        I = cls.sym2cat.assign(v=1).pivot(index="symbol", columns="goid", values='v')
        I = I.reindex(genes)
        I = I.fillna(0).astype(int)
        return I

    @classmethod
    def get_go_terms(cls):
        return pd.Series(nx.get_node_attributes(cls.G, name="name"), name="go_name")


def make_reduced_dataset():
    with Setup.download(Setup.URLS['meta']).now.open() as fd:
        log.info("Reading the full metadata.")
        df_meta = pd.read_csv(fd, sep=',', index_col=0)
        nsamples = len(df_meta)
        log.info(f"Original number of samples: {nsamples}.")

    if False:
        log.info("Creating a new reduced expression dataset.")
        with Setup.download(Setup.URLS['expr']).now.open() as fd:
            chunksz = 1024
            nchunks = (nsamples // chunksz) + 1

            from tcga.utils import seek_then_rewind
            with seek_then_rewind(fd):
                columns = pd.read_csv(fd, sep=',', nrows=0, index_col=0).columns

            columns = sorted(set(columns) & set(Setup.genes))

            objs = (
                df[columns].sample(axis=0, frac=Setup.frac_samples)
                for df in progressbar(pd.read_csv(fd, chunksize=chunksz, sep=',', index_col=0), max_value=nchunks)
            )

        # genes x samples
        df_expr = pd.concat(objs)
        filename = Setup.path_to_reduced / "data.csv.gz"
        df_expr.to_csv(filename, compression='gzip', sep='\t')
    else:
        log.info("Reading the existing reduced expression dataset.")
        df_expr = pd.read_table(Setup.path_to_reduced / "data.csv.gz", compression='gzip', index_col=0)

    # samples x meta
    log.info(f"Writing the reduced metadata ({len(df_expr)} samples).")
    df_meta = df_meta.loc[df_expr.index]
    filename = Setup.path_to_reduced / "meta.csv.gz"
    df_meta.to_csv(filename, compression='gzip', sep='\t')


def make_exploratory_plots():
    read = First(Setup.path_to_reduced.glob).then(unlist1).then(lambda f: pd.read_table(f, index_col=0))
    (df_meta, df_expr) = (read("meta*"), read("data*").T)

    I = Setup.get_incidence_matrix(df_expr.index)

    out_dir = mkdir(Path(__file__).with_suffix('') / "exploratory")

    def annotate(bars):
        # https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html
        for b in bars:
            px.a.annotate(
                f'{b.get_height()}',
                xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points", ha='center', va='bottom',
                fontsize=6,
            )

    def yscale_log(axis):
        axis.set_yscale('log')
        axis.set_ylim(0.91, 1.09 * (10 ** np.ceil(np.log10(max(px.a.get_ylim())))))
        axis.yaxis.set_major_formatter(lambda y, pos: f"{y:.0f}")

    with Plox() as px:
        gg = pd.Series(I.sum(axis=1)).value_counts()
        annotate(px.a.bar(gg.index, gg))
        px.a.set_xticks(gg.index)
        px.a.set_xlabel("... are in how many GO categories")
        px.a.set_ylabel("How many genes ...")
        yscale_log(px.a)
        px.f.savefig(out_dir / "go_by_gene.png")

    with Plox({rcParam.Xtick.labelsize: 5}) as px:
        cc = pd.Series(I.sum(axis=0)).value_counts().sort_index()
        rg = range(len(cc.index))
        annotate(px.a.bar(rg, cc))
        px.a.set_xticks(rg)
        px.a.set_xticklabels(cc.index)
        px.a.set_xlabel("... have how many genes")
        px.a.set_ylabel("How many GO categories ...")
        yscale_log(px.a)
        px.f.savefig(out_dir / "genes_by_go.png")


def main():
    import tensorflow as tf
    import tensorflow_addons as tfa

    # Random seed
    rs = 43

    # float64 doesn't work on model.fit for some reason
    TF_DTYPE = tf.float32

    #
    latent_n_per_goid = 1

    read = First(Setup.path_to_reduced.glob).then(unlist1).then(lambda f: pd.read_table(f, index_col=0))
    (df_meta, df_expr) = (read("meta*"), read("data*"))

    if df_expr.columns.equals(df_meta.index):
        df_expr = df_expr.T

    assert df_expr.index.equals(df_meta.index), "Each row should be a sample as in the metadata."

    log.warning(f"Removing samples with zero expression.")
    df_expr = df_expr[df_expr.any(axis=1)]
    log.info(f"Number of samples left: {len(df_expr)}")

    # Equilibrate training weight of each class
    # SAMPLE_WEIGHT_BY = 'cell_type_alias_label'
    # SAMPLE_WEIGHT_BY = 'subclass_label'
    SAMPLE_WEIGHT_BY = 'class_label'

    # log.warning(f"Normalizing sample-wise.")
    # assert all(df_expr.sum(axis=1))
    # df_expr = df_expr.div(df_expr.sum(axis=1), axis=0)

    # genes x goid incidence matrix
    I = Setup.get_incidence_matrix(df_expr.columns)
    assert all(pd.Series(I.columns).str.startswith("GO:"))

    sample_weight = (
        lambda s:
        # Make a series of weights ...
        pd.Series(name='w', index=s.index, data=(1 / s.value_counts(dropna=False)).reindex(s).tolist())
    )(
        # ... from this one
        df_meta[SAMPLE_WEIGHT_BY]
    )

    assert all(~sample_weight.isna())
    assert all(np.isclose(1, sample_weight.groupby(df_meta[SAMPLE_WEIGHT_BY]).sum()))

    batch_size = 16

    # print(tf.data.Dataset.from_tensor_slices(tf.cast(df_expr, dtype=TF_DTYPE)).take(1))
    # <TakeDataset shapes: (1033,), types: TF_DTYPE>

    ds_expr = tf.data.Dataset.from_generator(
        generator=df_expr.values.__iter__,  # iterate over rows
        output_signature=tf.TensorSpec(shape=(len(df_expr.columns),), dtype=TF_DTYPE),
    )

    # sample weight
    ds_splw = tf.data.Dataset.from_generator(
        generator=sample_weight.__iter__,
        output_signature=tf.TensorSpec(shape=(), dtype=TF_DTYPE)
    )

    # print(pd.Series(list(ds_smwt.as_numpy_iterator())).value_counts())

    ds = tf.data.Dataset.zip((ds_expr, ds_expr, ds_splw))
    ds = ds.batch(batch_size)
    ds = ds.shuffle(buffer_size=65, seed=43).prefetch(2)

    log.info("Dropping empty GO categories.")
    I = I[I.columns[I.any(axis=0)]]

    (n_genes, n_goids) = I.shape
    assert (n_genes == first(first(ds)).shape[1])

    with Peek(reporter=None) as peek:
        # feature extractors: goid-specific gene expression
        fex = [
            tf.SparseTensor(
                indices=peek(list(enumerate(first(np.where(I[goid]))))),
                values=peek(tf.constant([1] * np.count_nonzero(I[goid]), dtype=TF_DTYPE)),
                dense_shape=peek([np.count_nonzero(I[goid]), len(I.index)]),
            )
            for goid in I.columns
        ]

    # TODO:
    # deal w row-normalization; use log-trafo only?
    # what is the correct / statistically meaningful loss?
    # how to factor out the category size
    # add spurious categories
    # clustering performance
    # how to validate? artificial dataset?

    class Encoder(tf.keras.layers.Layer):
        def __init__(self, name="encoder"):
            super().__init__(name=name)

            trainable = True

            self.to_one = [
                self.add_weight(
                    name=f"to_one_{goid}",
                    shape=[f.shape[0], latent_n_per_goid],
                    initializer=tf.keras.initializers.RandomUniform(minval=0.01, maxval=0.1, seed=rs),
                    trainable=trainable,
                    dtype=TF_DTYPE,
                )
                for (goid, f) in zip(I.columns, fex)
            ]

        def call(self, x, training=False):
            # https://programming.vip/docs/day-6-tensorflow2-model-subclassing-api.html

            # log.info(x)

            x = tf.convert_to_tensor(x)

            norm1 = tf.reduce_sum(x, axis=1, keepdims=True)
            x = x / norm1

            fx = [
                tf.matmul(
                    tf.transpose(tf.sparse.sparse_dense_matmul(sp_a=f, b=x, adjoint_b=True)),
                    m
                )
                for (f, m) in zip(fex, self.to_one)
            ]

            #
            z = tf.concat(fx, axis=1)

            # https://tensorflow.google.com/api_docs/python/tf/debugging/Assert
            # this has no .mark_used()
            tf.debugging.Assert(z.shape[-1] == len(I.columns), data=[z])

            # z = tfa.layers.Sparsemax()(z)
            z = tf.keras.layers.Softmax()(z)

            output = {
                'go': z,
                'norm1': norm1,
            }

            return output

    class Decoder(tf.keras.layers.Layer):
        def __init__(self, name="decoder"):
            super().__init__(name=name)
            # input = tf.keras.layers.Input(shape=[n_genes])

            trainable = True

            self.un_one = [
                self.add_weight(
                    name=f"un_one_{goid}",
                    shape=[latent_n_per_goid, f.shape[0]],
                    initializer=tf.keras.initializers.RandomUniform(minval=0.01, maxval=0.1, seed=rs),
                    constraint=tf.keras.constraints.NonNeg(),
                    trainable=trainable,
                    dtype=TF_DTYPE,
                )
                for (goid, f) in zip(I.columns, fex)
            ]

        def call(self, latent, training=False):
            z = latent['go']
            norm1 = latent['norm1']

            tf.debugging.Assert(z.shape[-1] == len(I.columns), data=[z])
            zz = tf.split(z, axis=1, num_or_size_splits=(z.shape[1] // latent_n_per_goid))

            r = tf.transpose(tf.add_n([
                tf.sparse.sparse_dense_matmul(sp_a=f, b=tf.matmul(z, m), adjoint_a=True, adjoint_b=True)
                for (f, m, z) in zip(fex, self.un_one, zz)
            ]))

            r = (r / tf.reduce_sum(r, axis=1, keepdims=True)) * norm1

            return r

    inputs = tf.keras.Input(shape=[n_genes])
    encoder = Encoder()
    decoder = Decoder()
    model = tf.keras.Model(inputs=inputs, outputs=decoder(encoder(inputs)))

    try:
        model.get_layer("encoder").call(df_expr.iloc[0:3, :].astype('float32'))
        model.get_layer("encoder").call([first(ds_expr)])
    except:
        log.exception("Oops.")

    # Could be interesting:
    # MeanAbsolutePercentageError
    # MeanSquaredLogarithmicError
    # CosineSimilarity

    # model.run_eagerly = True

    opt = tf.keras.optimizers.Adam(learning_rate=1e-2)
    model.compile(optimizer=opt, loss=tf.keras.losses.MeanAbsolutePercentageError())

    # https://keras.io/guides/writing_your_own_callbacks/
    class MyCallback(tf.keras.callbacks.Callback):
        def on_epoch_end(self, epoch, logs=None):
            # log.info(f"End epoch {epoch} of training; got log keys: {list(logs.keys())}.")

            try:
                log.info(f"Writing Encoder(dataset)...")
                outfile = logpath / "latent.csv.gz"
                pd.DataFrame(
                    data=model.get_layer("encoder").call(tf.cast(df_expr, TF_DTYPE)).get('go').numpy(),
                    index=df_expr.index,
                    columns=I.columns,
                ).to_csv(
                    outfile, compression='gzip', sep='\t',
                )
                log.info(f"OK: {relpath(outfile)} .")
            except:
                log.exception("Oops.")

    callbacks = [
        tf.keras.callbacks.CSVLogger(logpath / "keras_log.csv"),
        tf.keras.callbacks.ModelCheckpoint(
            filepath=(logpath / "model"),
            save_weights_only=True,
            monitor='loss', mode='min',
            save_best_only=True,
        ),
        MyCallback(),
    ]

    history = model.fit(
        x=ds.shuffle(buffer_size=100).repeat(10),
        batch_size=None,
        #
        callbacks=callbacks,
        verbose=1,
        #
        shuffle=False,
        #
        # only for Tensors:
        validation_split=None,
        #
        # set if validation_data is provided as a tf.data:
        validation_steps=None,
        #
        initial_epoch=0,
        epochs=199,
        # steps_per_epoch=10,
        #
        # used for generator or keras.utils.Sequence input only:
        use_multiprocessing=True,
        workers=4,
    ).history

    with (logpath / "history.txt").open(mode='w') as fd:
        print(history, file=fd)


if __name__ == '__main__':
    # make_reduced_dataset()
    # make_exploratory_plots()
    main()

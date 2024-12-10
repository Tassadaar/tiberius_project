import sys, os
# ADD PATH TO LEARN MSA HERE (if not installed with pip)
learn_msa_path = '../../../install/learnMSA/'
sys.path.insert(0, learn_msa_path)

sys.path.append("../../bin")
import tensorflow as tf
from tensorflow.keras.optimizers import Adam
import tensorflow.keras as keras
from bin.eval_model_class import PredictionGTF
from bin.models import lstm_model, custom_cce_f1_loss
import matplotlib.pyplot as plt

batch_size = 20
seq_len = 9999
strand = '+'

emb=False
hmm_factor=1
trans_lstm=False

inp_data_dir = 'inp/'
out_dir = 'test_train/'
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

genome_path = f'{inp_data_dir}/genome.fa'
annot_path= f'{inp_data_dir}/annot.gtf'


pred_gtf = PredictionGTF(
    seq_len=seq_len,
    batch_size=batch_size,
    hmm=True,
    emb=False,
    num_hmm=1,
    hmm_factor=1,
    genome_path=genome_path,
    annot_path=annot_path,
    softmask=True, strand=strand,
)

# load input data x_seq
x_seq, y_seq, _ = pred_gtf.load_inp_data(
    strand=strand,
    chunk_coords=True, softmask=True, pad=False
)

# see lstm_model documentation for more arguments
config = {
    "num_epochs": 10,
    "stride": 0,
    "units": 100,
    "filter_size": 32,
    "numb_lstm": 2,
    "numb_conv": 3,
    "dropout_rate": 0.0,
    "hmm_dense": 32,
    "pool_size": 9,
    "lr": 1e-4,
    "batch_size": batch_size,
    "w_size": seq_len,
    'output_size': 15,
}

relevant_keys = ['units', 'filter_size', 'kernel_size',
                 'numb_conv', 'numb_lstm', 'dropout_rate',
                 'pool_size', 'stride',
                 'output_size', 'multi_loss']

relevant_args = {key: config[key] for key in relevant_keys if key in config}
model = lstm_model(**relevant_args)
adam = Adam(learning_rate=config["lr"])
f1loss = custom_cce_f1_loss(2, batch_size=config["batch_size"])


model.compile(loss=f1loss, optimizer=adam, metrics=['accuracy'])
model.summary()

history = model.fit(x=x_seq, y=y_seq,
          epochs=config["num_epochs"],
          batch_size=config["batch_size"])

# Save the trained model
model_save_path = os.path.join(out_dir, 'tiberius_lstm.h5')
model.save(model_save_path)
print(f"Model saved to {model_save_path}")

# Assuming you have 'history' from model.fit(...)
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('Training and Validation Loss')
plt.legend()

# Save the figure to a file instead of showing it
plt.savefig('training_loss_plot.png', dpi=300)  # You can adjust the dpi for resolution
plt.close()
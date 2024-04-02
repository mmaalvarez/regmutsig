# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# Autoencoder
# -----------
# - Adapted from Mischan Vali-Pour's Variational Autoencoder: https://github.com/lehner-lab/RDGVassociation/tree/main/somatic_component_extraction
# - ...which is in turn adapted/inspired from https://github.com/greenelab/tybalt/blob/master/tybalt_vae.ipynb
# - See also https://blog.keras.io/building-autoencoders-in-keras.html
#
# Main changes 
# ---
# - Now it's just a vanilla autoencoder, not "variational"
# - It uses a different permuted table (split into train and validation) in each epoch -- but if *n* epochs > *p* permuted tables (currently the case), it will cycle through them every *p* epochs

# +
import sys, os, argparse
import numpy as np
import pandas as pd
from glob import glob
# python ≥3.5 required
print(sys.version)

# tensorflow 1.15.5 recommended
import tensorflow as tf
print(tf.__version__)
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers
from tensorflow.keras.callbacks import Callback
from tensorflow.keras.utils import plot_model


# +
######## functions and classes

## to select optimizer
def get_optimizer(name, learning_rate):
    if name.lower() == 'adam':
        return optimizers.Adam(lr=learning_rate)
    elif name.lower() == 'sgd':
        return optimizers.SGD(lr=learning_rate)
    elif name.lower() == 'rmsprop':
        return optimizers.RMSprop(lr=learning_rate)
    else:
        raise ValueError('Optimizer name not recognized')
        
## custom callback to change the permuted training and validation data pair to use in each epoch          
class DataSwitchCallback(Callback):
    def __init__(self, data_dict):
        super(DataSwitchCallback, self).__init__()
        # initialize variables at first epoch only
        self.data_dict = data_dict # store train_validation_dfs_dict
        self.epoch_count = 0
    # when starting an epoch...
    def on_epoch_begin(self, epoch, logs=None):
        # ...find out which permuted data is to be used in this epoch
        data_index = self.epoch_count % len(self.data_dict)
        # ...select that permuted data (train and validation)
        self.model.train_data, self.model.val_data = self.data_dict[data_index]
        # ...print log message
        print(f"Using data pair #{data_index} for epoch {self.epoch_count + 1}/{epochs}")
        # ...and add another epoch to the count
        self.epoch_count += 1


# +
######## input arguments

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--dataset_real', 
                    default='autoencoder_input/original_scaled.tsv',
                    help='REAL dataset (i.e. no permutations), just [-1,1]-scaled, type directory + file name')
parser.add_argument('-p', '--dataset_permuted', 
                    default="autoencoder_input/permuted_coefficients_*__epoch_*.tsv",
                    help='PERMUTED training+validation samples (output from 1_generate_..._validating.R), type "DIRPATH/permuted_coefficients_*_epoch_*.tsv". IMPORTANT: it has to always go with double quotes in the python command so that the shell does not expand the wildcards, e.g. python -p "$dataset_permuted", or python -p "autoencoder_input/permuted_coefficients_*__epoch_*.tsv"')
parser.add_argument('-v', '--validation',
                    default='0.1',
                    help='random fraction of the dataset_permuted to keep aside as a validation set; required only for hyperparameter optimization')
parser.add_argument('-n', '--num_components',
                    default='5',
                    help='IMPORTANT to specify: latent space dimensionality (k, size, signatures)')
parser.add_argument('-d', '--depth',
                    default='1',
                    help='if depth=1 (default) there is only input, latent, and ouput layers; if depth=2, there will be 2 extra (hidden) layers: 1 between the input and the encoded layers, and 1 between the latent and the ouput layers')
parser.add_argument('-hd', '--hidden_dim',
                    default='2',
                    help='however many times the size of the encoded layer to be used as the size for the 2 extra hidden layers (only if depth==1)')
parser.add_argument('-b', '--batch_size',
                    default='200',
                    help='number of samples to include in each learning batch')
parser.add_argument('-e', '--epochs',         
                    default='200',
                    help='how many times to cycle -- every epoch a different set of PERMUTED training+validating samples is used')
parser.add_argument('-a', '--activation_function',         
                    default='relu',
                    help='activation function for all layers except the output (tanh hardcoded)')
parser.add_argument('-w', '--weight_initializer',         
                    default='glorot_uniform',
                    help='distribution from which initial weights are sampled')
parser.add_argument('-lf', '--loss_function',
                    default='mean_squared_error',
                    help='loss function')
parser.add_argument('-lr', '--learning_rate',
                    default='0.0005',
                    help='learning rate of the Adam optimizer')
parser.add_argument('-o', '--optimizer',
                    default='Adam',
                    help='optimizer algorithm')

# if interactive, pass values manually
if 'ipykernel' in sys.modules:
    args = parser.parse_args(['-r', 'autoencoder_input/original_scaled.tsv', \
                              '-p', "autoencoder_input/permuted_coefficients_*__epoch_*.tsv", \
                              '-v', '0.1', \
                              '-n', '5', \
                              '-d', '1', \
                              '-hd', '2', \
                              '-b', '200', \
                              '-e', '200', \
                              '-a', 'relu', \
                              '-w', 'glorot_uniform', \
                              '-lf', 'mean_squared_error', \
                              '-lr', '0.0005', \
                              '-o', 'Adam'])
else:
    args = parser.parse_args()

# +
######## set data paths and hyper-parameters
## removed the kappa and beta since they are specific to VAEs

dataset_real = args.dataset_real
dataset_permuted_list = glob(args.dataset_permuted) # this is a list of permuted tables
validation_set_fraction = float(args.validation)
latent_dim = int(args.num_components)
depth = int(args.depth)
hidden_dim = int(latent_dim)*int(args.hidden_dim) # set the 2 extra hidden layers (only added if depth==2) to be each a multiple/fraction of the size of the encoded layer
batch_size = int(args.batch_size)
epochs = int(args.epochs)
activation = args.activation_function
weight_init = args.weight_initializer
loss_function = args.loss_function
learning_rate = float(args.learning_rate)
optimizer_name = args.optimizer
optimizer = get_optimizer(optimizer_name, learning_rate)

# set random seed
seed = int(np.random.randint(low=0, high=10000, size=1))
np.random.seed(seed)

# +
######## upload and process input data

## upload [-1,1]-scaled ORIGINAL data (for final signature extraction)
real_df = pd.read_csv(dataset_real, sep='\t')
# store sample names column, renamed as "Sample"
sample_id = real_df.drop(real_df.columns[1:], axis=1).rename(columns={'sample_id': 'Sample'})
# store column names
column_names = real_df.drop(real_df.columns[0], axis=1).columns
# convert to numpy array (remove sample (first) column)
real_df = np.array(real_df.drop(real_df.columns[0], axis=1))
print(real_df.shape)

## upload [-1,1]-scaled full PERMUTED data (for training and validation)
# load all permuted tables, split each one into train and validation, and store in a dict; each pair will be used for a different epoch (or if n epochs > p permuted tables, just reuse them every p epochs) 
train_validation_dfs_dict = dict()
for i,dataset_permuted in enumerate(dataset_permuted_list):
    # load permuted table
    full_df = pd.read_csv(dataset_permuted, sep='\t')
    # split it into training and validation tables, using the specified fraction for validation 
    validation_df = full_df.sample(frac=validation_set_fraction, random_state=seed)
    train_df = full_df.drop(validation_df.index)
    # convert to numpy arrays (remove sample (first) and nIter (last) columns)
    train_df = np.array(train_df.drop(train_df.columns[0], axis=1).drop(train_df.columns[-1], axis=1))
    validation_df = np.array(validation_df.drop(validation_df.columns[0], axis=1).drop(validation_df.columns[-1], axis=1))
    # store
    train_validation_dfs_dict[i] = [train_df, validation_df]
    
# extract the 1st train/validation pair to know the number of neurons needed in input layer ("original_dim")
first_train_df,first_validation_df = train_validation_dfs_dict[0]
original_dim = first_train_df.shape[1]
print(first_train_df.shape)
print(first_validation_df.shape)

# +
#### ENCODER ####

# set number of neurons in input layer: i.e. 1 sample's coefficients vector, with 1 coefficient per predicting feature (DNA repair marks & SBS96)
coefficients_vector_input = Input(shape=(original_dim,))

# build hidden layers (at least the latent space, i.e. encoded layer)
if depth == 1: # default, only 1 encoded layer
    encoded = Dense(latent_dim, 
                    activation=activation, 
                    kernel_initializer=weight_init)(coefficients_vector_input)
elif depth == 2: # 1 hidden + the encoded layers
    hidden_dense = Dense(hidden_dim, 
                         activation=activation,
                         kernel_initializer=weight_init)(coefficients_vector_input)
    encoded = Dense(latent_dim, 
                    activation=activation, 
                    kernel_initializer=weight_init)(hidden_dense)

# +
#### DECODER ####

# build at least the output ("decoded") layer
if depth == 1: # default
    decoded = Dense(original_dim, 
                    # tanh hardcoded as activation function since our output values have to be between -1 and 1, like the input
                    activation='tanh', 
                    kernel_initializer=weight_init)(encoded)
elif depth == 2: # 1 hidden + 1 output layers
    hidden_decoded = Dense(hidden_dim, 
                           activation=activation, 
                           kernel_initializer=weight_init)(encoded)
    decoded = Dense(original_dim, 
                    activation='tanh', 
                    kernel_initializer=weight_init)(hidden_decoded)

# +
############### build autoencoder and encoder models ###############
## the 'autoencoder' and 'encoder' models are "intertwined through shared layer references. When you train the autoencoder, you're indirectly training the encoder because they share the same initial layers, and these layers have associated weights and biases. Since the encoder shares these layers with the autoencoder, any updates to the weights and biases in the autoencoder are reflected in the encoder. So after training, you can use the encoder separately to obtain the compressed (encoded) representations of data without reconstructing it, even though the training process involved both encoding and decoding"

## autoencoder model: map the input layer to its reconstruction (i.e. to the output layer)
autoencoder = Model(coefficients_vector_input, decoded)

# use optimizer to backpropagate based on the loss
autoencoder.compile(optimizer=optimizer, 
                    loss=loss_function) # mean squared error as loss function; 'binary_crossentropy' would yield higher losses

# Add a reference to the training and validation data in the autoencoder model
autoencoder.train_data = first_train_df
autoencoder.val_data = first_validation_df

# summary of autoencoder model
autoencoder.summary()


## encoder model: map the input layer to its encoded representation (i.e. to the latent space)
# this will be used after training
encoder = Model(coefficients_vector_input, encoded)

# +
############### train (fit) the 'autoencoder' model ##########################

# Create an instance of the callback with the data dictionary
data_switcher = DataSwitchCallback(train_validation_dfs_dict)

# Run the training
history = {'loss': [], 'val_loss': []}
for epoch in range(epochs):
    # use a different train_data + val_data pair at each iteration (i.e. epoch)
    epoch_it = autoencoder.fit(autoencoder.train_data, autoencoder.train_data,
                               shuffle = True,
                               epochs = 1, # run only one epoch at each loop iteration
                               batch_size = batch_size,
                               validation_data = (autoencoder.val_data, autoencoder.val_data),
                               callbacks = [data_switcher])
    # keep track of losses
    history['loss'].append(epoch_it.history['loss'][0])
    history['val_loss'].append(epoch_it.history['val_loss'][0])

# evaluate final loss
training_loss = autoencoder.evaluate(autoencoder.train_data, autoencoder.train_data)
validation_loss = autoencoder.evaluate(autoencoder.val_data, autoencoder.val_data)
print(f'Final training loss: {training_loss.round(2)}')
print(f'Final validation loss: {validation_loss.round(2)}')

# visualize training performance
history_df = pd.DataFrame(history)
ax = history_df.plot()
ax.set_xlabel('Epochs')
ax.set_ylabel('Autoencoder loss')
fig = ax.get_figure()

# +
##### now use the 'encoder' model (already trained via the autoencoder, see above) to encode the original (i.e. not permuted) -1,1 scaled coefficients matrix into the latent representation
## these k values should roughly correspond to each sample's k signature exposures from the NMF

encoded_real_df = encoder.predict_on_batch(real_df)

# rename columns ('signatures'), and convert to pandas df
encoded_real_df = pd.DataFrame(encoded_real_df, columns = range(1, latent_dim+1)).add_prefix('ae')

# append sample names column
encoded_real_df = pd.concat([sample_id, encoded_real_df], axis=1)

## check nodes activity to ensure that the model is learning a distribution of feature activations, and not zeroing out features
sum_node_activity = encoded_real_df.drop('Sample',axis=1).sum(axis=0).sort_values(ascending=False)
sum_node_activity_mean = sum_node_activity.mean()
print(sum_node_activity)

# +
############### save outputs ###############

output_folder_name = f'{str(validation_set_fraction)}-validation_' \
                     f'{str(latent_dim)}-components_' \
                     f'{str(depth)}-depth_' \
                     f'{str(hidden_dim)}-hiddendim_' \
                     f'{str(batch_size)}-batchsize_' \
                     f'{str(epochs)}-epochs_' \
                     f'{activation}_' \
                     f'{weight_init}_' \
                     f'{loss_function}_' \
                     f'{str(learning_rate)}-learningrate_' \
                     f'{optimizer_name}_' \
                     f'{str(round(sum_node_activity_mean,2))}-meansumactivity_' \
                     f'{str(seed)}-seed/'

## if interactive
if 'ipykernel' in sys.modules:
    if not os.path.exists('autoencoder_output'):
        os.mkdir('autoencoder_output')
    output_folder_name = 'autoencoder_output/' + output_folder_name
    os.mkdir(output_folder_name)
else:
    ## not interactive (nextflow handles the autoencoder_output folder creation)
    os.mkdir(output_folder_name)
    
# training performance plot
fig.savefig(output_folder_name + 'loss_history.jpg', dpi=600)    

# viz network
plot_model(autoencoder,
           to_file=output_folder_name + 'model_viz.jpg',
           show_layer_names=False,
           show_shapes=True,
           rankdir='TB', #"LR" horizontal
           dpi=600)

# encoder model
encoder.save(output_folder_name + 'encoder_model.tf')

# encoded layer ("signature exposures")
encoded_real_df.to_csv(output_folder_name + 'encoded_layer.tsv', sep='\t', index= False)

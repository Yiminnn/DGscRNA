
import sys
import scanpy
import os
import torch
import numpy as np
import torchmetrics
import numpy as np
import loompy
import pandas as pd
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data_utils
from torch.autograd import Variable

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)




def get_stats(preds, labels, num_possible_labels):
    
    # get all stats
    all_stats = {}
    
    # accuracy
    acc = torchmetrics.Accuracy(task = 'multiclass', num_classes = num_possible_labels)
    all_stats['acc'] = acc(preds, labels).item()
    
    # F1Score
    f1 = torchmetrics.F1Score(task = 'multiclass', num_classes = num_possible_labels)
    all_stats['f1'] = f1(preds, labels).item()
    
    # PR
    all_stats['p'] = torchmetrics.Precision(task = 'multiclass', num_classes = num_possible_labels)(preds, labels).item()
    all_stats['r'] = torchmetrics.Recall(task = 'multiclass', num_classes = num_possible_labels)(preds, labels).item()
    
    return(all_stats)

# DGCyTOF DeepModel
class DeepModel(nn.Module):
    def __init__(self, input_shape, class_shape):
        super().__init__()
        self.activation = nn.LeakyReLU()
        self.output_activation = nn.Softmax()
        self.fc1 = nn.Linear(input_shape, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, class_shape)

    def forward(self, x):
        x = self.activation(self.fc1(x))
        x = self.activation(self.fc2(x))
        x = self.output_activation(self.fc3(x))
        
        return x

def get_num_correct(preds, labels):
    return preds.argmax(dim=1).eq(labels).sum().item()

def run_dgscrna(folder_loc, data_name, study_sets, EPOCHS = 10):
    # read the data
    data_w_batch_corr_annotated_loc = os.path.join(folder_loc, data_name)
    data_w_batch_corr_annotated = scanpy.read_loom(data_w_batch_corr_annotated_loc, sparse = False)

    # get annotations
    orig_annotations = data_w_batch_corr_annotated.obs.copy()

    # get analysis columns (has + in the name)
    analyze_columns = [i for i in orig_annotations.columns if ('_none' in i) or ('_0.5' in i) or ('_mean' in i)] #('clusters' in i) and 

    orig_annotations.index = [i.split('_')[0] for i in orig_annotations.index]
    orig_annotations.index = orig_annotations.index.values + '_' + np.arange(start = 0, stop = orig_annotations.shape[0]).astype(str)

    orig_annotations['cell_id'] = np.arange(0, orig_annotations.shape[0])

    # load results if present, if not continue
    annotation_file = 'annotated_' + data_name + '.csv'
    if os.path.isfile(annotation_file):
        final_res = pd.read_csv(annotation_file, index_col = 0)
    else:
        final_res = orig_annotations.copy()
        final_res.to_csv(annotation_file)

    # proceed to run DGscRNA on all columns
    print('Entering Loop')
    for different_annotation in analyze_columns:
        
        print('Working on ' + different_annotation)
        
        # get cell types
        final_cell_types = orig_annotations[different_annotation]
        print(final_cell_types.value_counts())

        if different_annotation + '_DGscRNA' in final_res.columns:
            continue

        if (final_cell_types == 'Undecided').sum() == final_cell_types.shape[0]:
            print('All undecided/unknown')
            final_res[different_annotation + '_DGscRNA'] = ['Undecided'] * final_cell_types.shape[0]
            final_res.to_csv(annotation_file)
            continue
            
        if (final_cell_types == 'Undecided').sum() == 0:
            print('No undecided/unknown')
            final_res[different_annotation + '_DGscRNA'] = final_cell_types
            final_res.to_csv(annotation_file)
            continue

        agreed_cells = final_cell_types.loc[final_cell_types != 'Undecided'].index
        disagreed_cells = final_cell_types.loc[final_cell_types == 'Undecided'].index

        print('Num agreed: ' + str(len(agreed_cells)))
        print('Num disagreed: ' + str(len(disagreed_cells)))

        all_cell_types = set(orig_annotations[different_annotation])

        bag_of_classes = {}
        rev_bag_of_classes = {}
        j = 0
        for i in sorted(list(all_cell_types)):
            bag_of_classes[i] = j
            rev_bag_of_classes[j] = i
            j += 1

        agreed_labels = [bag_of_classes[i] for i in orig_annotations.loc[agreed_cells, different_annotation]]
        non_agreed_labels = [bag_of_classes[i] for i in orig_annotations.loc[disagreed_cells, different_annotation]]

        agreed_cells_idx = orig_annotations.loc[orig_annotations[different_annotation] != 'Undecided', 'cell_id']
        disagreed_cells_idx = orig_annotations.loc[orig_annotations[different_annotation] == 'Undecided', 'cell_id']

        ###### prepare data

        y_train = torch.tensor(agreed_labels, dtype=torch.long)
        X_train = torch.tensor(data_w_batch_corr_annotated.X[agreed_cells_idx, ], dtype=torch.float) 

        y_test = torch.tensor(non_agreed_labels, dtype=torch.long)
        X_test = torch.tensor(data_w_batch_corr_annotated.X[disagreed_cells_idx, ], dtype=torch.float) 


        train_tensor = data_utils.TensorDataset(X_train, y_train) 
        test_tensor = data_utils.TensorDataset(X_test, y_test) 

        train_size = round(len(agreed_cells) * 0.90)
        val_size = len(agreed_cells) - train_size

        train_tensor, val_tensor = torch.utils.data.random_split(train_tensor, [train_size, val_size], generator=torch.Generator().manual_seed(42))

        params = {"batch_size": 256,
            "shuffle": False,
            "num_workers": 8}

        train_data = data_utils.DataLoader(train_tensor, **params) 
        test_data = data_utils.DataLoader(test_tensor, **params) 
        val_data = data_utils.DataLoader(val_tensor, **params) 

        del train_tensor
        del test_tensor
        del val_tensor

        print('----------------------------------------------')

        model = DeepModel(input_shape = 2000, class_shape = len(bag_of_classes))

        lr = 1e-3
        batch_size = params["batch_size"]
        criterion = nn.CrossEntropyLoss()
        optimizer = torch.optim.Adamax(model.parameters(), lr=lr)

        # train the data

        for epoch in range(EPOCHS):  # loop over the dataset multiple times

            total_loss = 0

            all_preds = []
            all_labels = []

            for batch in train_data: # Get Batch
                samples, labels = batch 

                preds = model(samples) # Pass Batch
                loss = criterion(preds, labels) # Calculate Loss

                optimizer.zero_grad() # Zero Gradients      
                loss.backward() # Calculate Gradients
                optimizer.step() # Update Weights

                total_loss += loss.item() * batch_size
                preds = preds.argmax(dim=1)

                all_preds.append(preds)
                all_labels.append(labels)

            all_preds = torch.cat(all_preds)
            all_labels = torch.cat(all_labels)
            all_stats = get_stats(all_preds, all_labels, num_possible_labels = len(bag_of_classes))

            print("epoch", epoch, "| ACC:", all_stats['acc'],  "| F1:", all_stats['f1'], "| Precision:", all_stats['p'], "| Recall:", all_stats['r'], "loss:", total_loss)

        print ('-'*80)   
        # tb.close()

        print('----------------------------------------------')

        # do validation

        classes = [i for i in bag_of_classes]
        class_correct = list(0. for i in range(len(bag_of_classes)))
        class_total = list(0. for i in range(len(bag_of_classes)))

        val_correct = 0
        val_total = 0

        for data in val_data:
            val_samples, val_labels = data
            val_outputs = model(Variable(val_samples))
            _, val_predicted = torch.max(val_outputs.data, 1)   # Find the class index with the maximum value.
            c = (val_predicted == val_labels).squeeze()
            for i in range(val_labels.shape[0]):
                label = val_labels[i]
                class_correct[label] += c[i].item()
                class_total[label] += 1

            val_total += val_labels.size(0)
            val_correct += (val_predicted == val_labels).sum()

        print("Accuracy:", round(100 *val_correct.item() / val_total, 4))
        print('-'*100)
        for i in range(len(bag_of_classes)):
            if class_total[i] == 0:
                print('Class >||{}||< not in set'.format(classes[i]))
            else:
                print('Accuracy of {} : {}'.format (
                classes[i], round(100 * class_correct[i] / class_total[i], 3)))

        # get stats
        all_val_preds = []
        all_val_preds_probs = []
        all_val_labels = []

        for batch in val_data: # Get Batch
            samples, labels = batch 

            preds = model(samples) # Pass Batch

            all_val_preds_probs.append(preds)

            preds = preds.argmax(dim=1)

            all_val_preds.append(preds)
            all_val_labels.append(labels)

        all_val_preds = torch.cat(all_val_preds)
        all_val_labels = torch.cat(all_val_labels)
        all_val_preds_probs = torch.cat(all_val_preds_probs)

        all_val_stats = get_stats(all_val_preds, all_val_labels, num_possible_labels = len(bag_of_classes))

        print("ACC:", all_val_stats['acc'],  "| F1:", all_val_stats['f1'], "| Precision:", all_val_stats['p'], "| Recall:", all_val_stats['r'])

        print('Checking for non available classes in validation: (prediction)')
        print([i for i in bag_of_classes if bag_of_classes[i] not in torch.unique(all_val_preds).tolist()])

        # apply on the unknowns
        all_test_preds = []
        all_test_preds_probs = []
        all_test_labels = []
        all_test_samples = []

        for batch in test_data: # Get Batch
            samples, labels = batch 

            preds = model(samples) # Pass Batch

            all_test_preds_probs.append(preds)

            preds = preds.argmax(dim=1)

            all_test_preds.append(preds)
            all_test_labels.append(labels)
            all_test_samples.append(samples)

        all_test_preds = torch.cat(all_test_preds)
        all_test_labels = torch.cat(all_test_labels)
        all_test_preds_probs = torch.cat(all_test_preds_probs)
        all_test_samples = torch.cat(all_test_samples)

        all_test_stats = get_stats(all_test_preds, all_test_labels, num_possible_labels = len(bag_of_classes))

        print("ACC:", all_test_stats['acc'],  "| F1:", all_test_stats['f1'], "| Precision:", all_test_stats['p'], "| Recall:", all_test_stats['r'])

        # get probabilities from valdiation
        probabilities = np.max(all_val_preds_probs.data.numpy(), axis=1)
        print('Minimum probability')
        print(min(probabilities))

        probabilities = np.max(all_test_preds_probs.data.numpy(), axis=1)
        tem = [round(i,4) for i in probabilities]
        correct_index = np.array([i for i,v in enumerate(tem) if v >= 0.90])
        incorrect_index = np.array([i for i,v in enumerate(tem) if v < 0.90])

        # print('correct index:')
        # print(len(correct_index))

        # print('incorrect index:')
        # print(len(incorrect_index))

        # Correct labeled ones are correct
        final_cell_types.iloc[disagreed_cells_idx.iloc[correct_index]] = [rev_bag_of_classes[i] for i in all_test_preds[correct_index].numpy()]
        final_cell_types.iloc[disagreed_cells_idx.iloc[incorrect_index]] = 'Unknown'

        # print('Final calculated types')
        # print(final_cell_types.value_counts())

        final_cell_types.replace('0', 'No_Annotation', inplace = True)

        # get the general 
        if study_sets:
            tissue_general = []
            for ct in final_cell_types:
                entry = str(ct).split('+')
                if 'HPA' in str(ct):# == 'HPA':
                    tissue_general.append(entry[0])
                else:
                    if entry[0] == 'Undecided' or entry[0] == 'Unknown' or entry[0] == 'No_Annotation':
                        tissue_general.append(entry[0])
                    elif entry[0] == '0':
                        tissue_general.append('No_Annotation')
                    elif (entry[0] == 'CellMarker_normal') or (entry[0] == 'CellMarker_cancer'):
                        tissue_general.append('+'.join(entry[3:len(entry)]))
                    elif entry[0] == 'NCOMMREFF':
                        tissue_general.append(entry[1])
                    else:
                        tissue_general.append('+'.join(entry[2:len(entry)]))

            final_res[different_annotation + '_DGscRNA_General'] = tissue_general

        final_res[different_annotation + '_DGscRNA'] = final_cell_types
        
        
        final_res.to_csv(annotation_file)

    

    

if __name__ == "__main__":
    # location where the data is
    folder_loc = str(sys.argv[1])
    # data file name
    data_name = str(sys.argv[2])
    # if custom sets or study sets have been used
    study_sets = str(sys.argv[3]) == 'study_sets'
    print('Start DGscRNA')
    run_dgscrna(folder_loc, data_name, study_sets)

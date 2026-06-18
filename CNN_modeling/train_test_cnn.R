train_test_cnn <- function(testing_label_vec, image_input_directory = getwd(),
                           architecture = 'VGG19',
                           model_output_filename = NULL,
                           validation_proportion = 0.1, n_epochs = 75,
                           fine_tune_epochs = NULL) {
  # Trains a convolutional neural network (CNN) based on the labels assigned, 
  # with labeled_chunks, to sub-areas of raster_layer.
  #
  # Args:
  #   testing_label_vec: A vector of class integer that provides the labels
  #     of the testing images. It is intended that this is the vector returned
  #     by the generate_images function.
  #   image_input_directory: A string identifying the directory to which the
  #     training and testing images have been written. It must have
  #     subdirectories called 'training_dir' and 'testing_dir,' and
  #     'training_dir' must, in turn, have two subdirectories, 'class_neg' and
  #     'class_pos.' It is intended that image_input_directory
  #     is the same directory specified by image_output_directory when running
  #     the generate_images function.
  #   architecture: A string indicating the CNN model architecture to be 
  #     implemented. Must be one of 'VGG16,' 'VGG19,' 'ResNet50,' 'ResNet152,'
  #     'InceptionV3,' or 'Xception.' 'VGG19' is default.
  #   model_output_filename: The desired file name for the trained model.
  #   validation_proportion: The proportion of the training data that will
  #     be used for validation while the model is being trained (this is not
  #     the same as testing data). Default is 0.1.
  #   n_epochs: The number of epochs for training the model. Default is 75.
  #   fine_tune_epochs: If not NULL, the model will be fine tuned for this
  #     number of epochs. NULL by default.
  #   
  # Returns:
  #   A 5 x 3 matrix containing measures of the model's performance on the
  #   testing data: precision, recall, F1 score, and area under the 
  #   receiver operating characteristic (ROC) curve (AUC); the trained model
  #   itself, though not returned by the function, is saved
  #   to model_output_filename.
  
  # Make sure required packages are installed.
  png_installed <- require(png)
  magrittr_installed <- require(magrittr)
  keras_installed <- require(keras)
  tfhub_installed <- require(tfhub)
  pins_installed <- require(pins)
  tfd_installed <- require(tfdatasets)
  roc_installed <- require(pROC)
  stopifnot(png_installed & magrittr_installed & keras_installed &
              tfhub_installed & pins_installed & tfd_installed & roc_installed)
  
  # Make sure valid name is provided to identify desired model architecture.
  possible_architectures <- c('VGG16', 'VGG19', 'ResNet50',
                              'ResNet152', 'InceptionV3', 'Xception')
  if (architecture %in% possible_architectures == FALSE) {
    stop('\'architecture\' must be \'VGG16\', \'VGG19\', or \'ResNet152\'')
  }
  
  # Make sure an output file name for the tf model has been specified.
  if (is.null(model_output_filename)) stop('Need output model file name.')
  
  # Identify the directories that hold the image files that will be used to
  # train and test the model.
  training_dir <- paste0(image_input_directory, '/training_dir')
  testing_dir <- paste0(image_input_directory, '/testing_dir')
  train_neg_dir <- paste0(image_input_directory, '/training_dir/class_neg')
  train_pos_dir <- paste0(image_input_directory, '/training_dir/class_pos')
 
  print('Start training CNN.')
  
  # Get image size.
  image_size <- dim(
    readPNG(
      paste0(
        image_input_directory, '/testing_dir/', dir(
          paste0(image_input_directory, '/testing_dir'))[1])))[1]
  
  # The following two functions read in the png image data in a manner that's
  # convenient for all we're about to do next. The documentation for the function
  # (image_datasets_from_directory) explains what's happening here nicely.
  training_ds <- image_dataset_from_directory(
    training_dir,
    labels = 'inferred',
    color_mode = 'rgb',
    image_size = c(image_size, image_size),
    seed = 2113,
    validation_split = validation_proportion,
    subset = 'training'
  )
  
  val_ds <- image_dataset_from_directory(
    training_dir,
    labels = 'inferred',
    color_mode = 'rgb',
    image_size = c(image_size, image_size),
    seed = 2113,
    validation_split = validation_proportion,
    subset = 'validation'
  )
  
  # Define a data augmentation function.
  data_augment <- keras_model_sequential() %>%
    layer_random_flip('horizontal') %>%
    layer_random_rotation(0.1)
  
  # Identify the manner in which the image data needs to be processed before
  # being ingested into the model. Most of the architecture options require
  # the 'caffe' preprocessing mode, though some will need 'tf.'
  preprocess_mode <- 'caffe'
  
  # The following if-else statements identify and load in the specified
  # CNN architecture, as well as the weights it's learned on
  # the IMAGENET data set. When loading in this model, though, the top layer is
  # excluded, which leaves a place for us to put another model on top.
  if (architecture == 'VGG16') {
    conv_base <- application_vgg16(
      weights = 'imagenet',
      include_top = FALSE,
      input_shape = c(image_size, image_size, 3)
    )
  } else if (architecture == 'VGG19') {
    conv_base <- application_vgg19(
      weights = 'imagenet',
      include_top = FALSE,
      input_shape = c(image_size, image_size, 3)
    )
  } else if (architecture == 'ResNet50') {
    conv_base <- application_resnet50(
      weights = 'imagenet',
      include_top = FALSE,
      input_shape = c(image_size, image_size, 3)
    )
  } else if (architecture == 'ResNet152') {
    conv_base <- application_resnet152(
      weights = 'imagenet',
      include_top = FALSE,
      input_shape = c(image_size, image_size, 3)
    )
  } else if (architecture == 'InceptionV3') {
    conv_base <- application_inception_v3(
      weights = 'imagenet',
      include_top = FALSE,
      input_shape = c(image_size, image_size, 3)
    )
    preprocess_mode <- 'tf'
  } else if (architecture == 'Xception') {
    conv_base <- application_xception(
      weights = 'imagenet',
      include_top = FALSE,
      input_shape = c(image_size, image_size, 3)
    )
    preprocess_mode <- 'tf'
  }
  
  conv_base$trainable <- FALSE  # Freeze the base model -- that is, prevent its
  # weights from being updated when, below, we run through it the training data
  # specific to our application. When we do that, we want to train only the
  # weights for the new model we're training (at least, for the moment).
  
  # Now we put a new model on top of the one that's already trained.
  inputs <- layer_input(c(image_size, image_size, 3))
  outputs <- inputs %>%
    data_augment() %>%
    {if (preprocess_mode == 'caffe') imagenet_preprocess_input(., mode = 'caffe') else imagenet_preprocess_input(., mode = 'tf')} %>%
    conv_base(training = FALSE) %>%  # This is not the same as freezing the 
    # existing model, which we did by setting conv_base$trainable to FALSE. 
    # This, instead, will keep the batch normalization layers in the base model
    # from updating, even when that model in general is unfrozen, which it will
    # be below.
    layer_global_average_pooling_2d() %>%
    layer_dropout(.2) %>%
    layer_dense(units = 1, activation = 'sigmoid')
  new_model <- keras_model(inputs, outputs)
  
  # Compile the model.
  new_model %>% compile(
    optimizer = optimizer_adam(),
    loss = loss_binary_crossentropy(),
    metrics = metric_binary_accuracy()
  )
  
  # This is where we train the model on our own (new) data.
  new_model %>% fit(
    training_ds,
    epochs = n_epochs,
    validation_data = val_ds
  )
  
  # Fine-tune the model if fine_tune_epochs is not NULL.
  if(!is.null(fine_tune_epochs)) {
    conv_base$trainable <- TRUE
    new_model %>% compile(
      optimizer = optimizer_adam(1e-3),  # Note how small the learning
      # rate we're setting here is.
      loss = loss_binary_crossentropy(),
      metrics = metric_binary_accuracy()
    )
    
    new_model %>% fit(
      training_ds,
      epochs = fine_tune_epochs,
      validation_data = val_ds
    )
  }
  save_model_tf(new_model, filepath = model_output_filename)
  print('Trained model saved.')
  
  print('Run predictions on testing data.')
  test_ds <- image_dataset_from_directory(
    testing_dir,
    labels = testing_label_vec,
    color_mode = 'rgb',
    image_size = c(image_size, image_size),
    shuffle = FALSE
  )
  
  # Calculate the model's loss and binary accuracy.
  test_results <- new_model %>% evaluate(test_ds)
  print('Loss and binary accuracy:')
  print(test_results)
  
  print('Generate other measures of model performance.')
  test_output_probs <- new_model %>% predict(test_ds)
  test_output_labels <- test_output_probs %>% `>`(0.5) %>% k_cast('int32')
  test_output_labels <- as.numeric(test_output_labels[, 1])
  
  results_df <- data.frame(predicted = test_output_labels,
                           true = testing_label_vec, 
                           result = vector(length = length(test_output_labels)))
  for (i in 1:nrow(results_df)) {
    if (results_df[i, 1] == 1) {
      if (results_df[i, 2] == 1) {
        results_df[i, 3] <- 'TP'
      } else {
        results_df[i, 3] <- 'FP'
      }
    } else {
      if (results_df[i, 2] == 0) {
        results_df[i, 3] <- 'TN'
      } else {
        results_df[i, 3] <- 'FN'
      }
    }
  }
  tps <- nrow(results_df[which(results_df$result == 'TP'), ])
  fps <- nrow(results_df[which(results_df$result == 'FP'), ])
  tns <- nrow(results_df[which(results_df$result == 'TN'), ])
  fns <- nrow(results_df[which(results_df$result == 'FN'), ])
  
  precision <- tps / (tps + fps)
  recall <- tps / (tps + fns)
  f1 <- 2 * (precision * recall) / (precision + recall)
  print(paste0('Precision: ', precision))
  print(paste0('Recall: ', recall))
  print(paste0('F1 score: ', f1))
  
  print('Generating ROC object.')
  test_roc <- roc(testing_label_vec, as.numeric(test_output_probs))
  print(test_roc)
  
  result_df <- data.frame(matrix(data = c('TP', 'FP', 'TN', 'FN', 'Bin_acc', 'Precision', 'Recall', 'F1', 'AUC',
                                   tps, fps, tns, fns, test_results[2], precision, recall, f1, auc(test_roc)),
                          ncol = 2))
  colnames(result_df) <- c('Measure', 'Score')
  result_df$Score <- as.numeric(result_df$Score)
  return(result_df)
}


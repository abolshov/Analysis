import tensorflow as tf
from NetUtils import TrainableSiLU, CombinedLoss


path_to_pretrained = "model_hh_dl_all_aug_trainable_silu_dynamic_comb_loss.keras"
custom_objects = {"TrainableSiLU": TrainableSiLU,
                  "CombinedLoss": CombinedLoss}
model = tf.keras.models.load_model(path_to_pretrained,
                                   custom_objects=custom_objects)
print(model.summary())
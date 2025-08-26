import tensorflow as tf
import tf2onnx
import onnx
# from onnx_tf.backend import prepare

from ModelUtils import LoadModel, ConvertToONNX

def main():
    model, _ = LoadModel('model_cfg.yaml')
    # print(model.name)
    print(model.summary())

    input_signature = [tf.TensorSpec([None, 157], tf.float32, name='input')]
    onnx_model = ConvertToONNX(model, input_signature, 'conversion_test')
    # onnx_model, _ = tf2onnx.convert.from_keras(model,
    #                                            input_signature=input_signature,
    #                                            opset=13)

    # onnx.save(onnx_model, 'conversion_test/model.onnx')

    onnx_model = onnx.load('conversion_test/model.onnx')
    onnx.checker.check_model(onnx_model)
    # print(onnx.helper.printable_graph(onnx_model.graph))

    # tf_rep = prepare(onnx_model)
    # tf_rep.export_graph('tf_model')

    # tf_model = tf.saved_model.load('tf_model')
    # infer = model.signatures['seving_default']

if __name__ == '__main__':
    main()
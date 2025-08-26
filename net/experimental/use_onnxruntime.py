import onnxruntime as ort
import numpy as np

def main():
    session = ort.InferenceSession('conversion_test/model.onnx')

    input_name = session.get_inputs()[0].name
    output_names = [out.name for out in session.get_outputs()]

    print(f'input_name: {input_name}')
    print(f'output_names: {output_names}')

    x = np.random.rand(4, 157).astype(np.float32)

    outputs = session.run(output_names, {input_name: x})

    for i, out in enumerate(outputs):
        print(f'output {i} shape: {out.shape}')

if __name__ == '__main__':
    main()

import h5py
import tensorflow as tf


def _floats_feature(value):
    return tf.train.Feature(float_list=tf.train.FloatList(value=value))

def _bytes_feature(value):
    """Returns a bytes_list from a string / byte."""
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))
    

def get_feature(h5_file, index):
    return {
        'features': _bytes_feature(h5_file['features'][index].flatten().tobytes()),
        'meta': _floats_feature(h5_file['meta'][index].flatten()),
        'labels_k': _bytes_feature(h5_file['labels_k'][index].flatten().tobytes()),
        'labels_r': _bytes_feature(h5_file['labels_r'][index].flatten().tobytes()),
        'probe_r': _bytes_feature(h5_file['probe_r'][index].flatten().tobytes()),
    }

def h5_to_tfrecord_converter(input_files, output_file_path):
    with tf.io.TFRecordWriter(output_file_path) as writer:
        for idx, (_input, _file_name) in enumerate(input_files):
            print('\n\ton job %d/%d, %s' % (idx, len(files), _file_name), end='')
            h5_file = h5py.File(_input)
            
            num_of_items = h5_file['meta'][:].shape[0]
            
            for index in range(num_of_items):
                example = tf.train.Example(
                features=tf.train.Features(
                    feature = get_feature(h5_file,index)
                ))
                writer.write(example.SerializeToString())
                print('\r{:.1%}'.format((index+1)/num_of_items), end='')
    

if __name__ == "__main__":

    import os

    INPUT_PATH = '/media/thomas/SSD/Data_64/Validation/'
    OUTPUT_PATH = ''
    RECURSIVE = True

    if INPUT_PATH.endswith('.h5'):
        if OUTPUT_PATH == '':
            OUTPUT_PATH = INPUT_PATH[:-3]
        print('Start converting...\t')
        h5_to_tfrecord_converter([INPUT_PATH], os.path.abspath(OUTPUT_PATH) + '.tfrecord')

    elif RECURSIVE:
        files = []
        if OUTPUT_PATH == '':
            OUTPUT_PATH = INPUT_PATH + 'data.tfrecord'
        for _file in os.listdir(INPUT_PATH):
            if _file.endswith('.h5'):
                files.append((os.path.join(INPUT_PATH, _file), _file[:-3] ))
        print(len(files), 'of HDF5 file detected.')
        h5_to_tfrecord_converter(files, OUTPUT_PATH)

    else:
        pass

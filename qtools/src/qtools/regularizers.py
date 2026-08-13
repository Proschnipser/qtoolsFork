import tensorflow as tf
from tensorflow.keras import regularizers

class WeightEntropyRegularizer(regularizers.Regularizer):
    def __init__(self, mean_vector, strength=5e-1, epsilon=1e-8, temperature=1):
        """
        mean_vector: 1D array/tensor, shape (input_dim,)
                     mean of each input position across all training vectors
        """
        self.mean_vector = tf.constant(mean_vector, dtype=tf.float32)
        self.strength = strength # equivalent to lambda
        self.epsilon = epsilon # to avoid log(0)
        self.temperature = temperature

    def __call__(self, kernel):
        # kernel shape: (input_dim, num_outputs) -> e.g. (1891, 5)
        # normalize each output neuron's weight column by the mean input vector
        normalized = kernel * (self.mean_vector[:, None])

        # softmax over input positions (axis=0) -> distribution per output neuron
        probs = tf.nn.softmax(normalized/self.temperature, axis=0)

        # Shannon entropy per output neuron (per column)
        entropy_per_output = -tf.reduce_sum(
            probs * tf.math.log(probs + self.epsilon), axis=0
        )

        mean_e = tf.reduce_mean(entropy_per_output)
        var_e  = tf.math.reduce_variance(entropy_per_output)

        return self.strength * (mean_e + mean_e*var_e)
        #return self.strength * tf.reduce_mean(entropy_per_output)

    def get_config(self):
        return {
            "mean_vector": self.mean_vector.numpy().tolist(),
            "strength": self.strength,
            "epsilon": self.epsilon,
            "temperature": self.temperature,
        }
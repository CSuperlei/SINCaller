from keras.preprocessing.text import one_hot
from keras.preprocessing.sequence import pad_sequences
from keras.models import Sequential
from keras.layers import Dense, LSTM, Dropout, Bidirectional
from keras.layers import Flatten
from keras.layers.embeddings import Embedding
from keras.utils.np_utils import to_categorical


# define documents
docs = ['a, t, t',
		'a, a, a',
		'g, t, t, t, t, t',
		'c, c, c, c, c, c',
		'a, g, g',
		't, t, t',
		't, c, c',
		'c, c, c',
		'g, c, c, c',
		'a, a, a, a',
		'g, g, g, t',
		't, t, t, t',
		'a, g, g, g, g, g, g, g, g, g, g, g, g, g, g',
		'g, g, g, g, g, g, g, g, g, g, g, g, g, g, g',
		'g a a a',
		'c c c c',
        'a g g g g g f f f g g g g g g f',
        'g g g g g g a a a g g g g g g g g',
        'a t c g',
        't t t t t t t t t t c a',
        'a a a a a t c',
        'c c a t c c c',
        'g a t c',
        'a a c t'
        ]
# define class labels
labels = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 2, 2, 2, 2, 2, 2, 2]
ohl=to_categorical(labels)
print(ohl)
# integer encode the documents
vocab_size = 20
encoded_docs = [one_hot(d, vocab_size) for d in docs]
print(encoded_docs)
# pad documents to a max length of 4 words
max_length = 15
padded_docs = pad_sequences(encoded_docs, maxlen=max_length, padding='post')
print(padded_docs)
print('padded_docs shape', padded_docs.shape)
# define the model
model = Sequential()
model.add(Embedding(vocab_size, 8, input_length=max_length))

# model.add(LSTM(128, return_sequences=True))
# model.add(Dropout(0.1))
# model.add(LSTM(256, return_sequences=False))
model.add(Bidirectional(LSTM(128, return_sequences=True)))
model.add(Bidirectional(LSTM(64, return_sequences=False)))
model.add(Dense(32))
model.add(Dense(3, activation='softmax'))
# compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
# summarize the model
print(model.summary())
# fit the model
model.fit(padded_docs, ohl, epochs=200, verbose=1)
# evaluate the model
loss, accuracy = model.evaluate(padded_docs, ohl, verbose=1)
print('Accuracy: %f' % (accuracy*100))
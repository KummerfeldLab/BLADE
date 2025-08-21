import pickle

def load_data():
    with open('data/sample.pkl', 'rb') as file:
        sample = pickle.load(file)
    return sample 
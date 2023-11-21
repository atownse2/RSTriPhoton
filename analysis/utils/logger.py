class Logger:
    def __init__(self, max_verbosity=0):
        self.max_verbosity = max_verbosity

    def set_max_verbosity(self, verbosity):
        self.max_verbosity = verbosity

    def log(self, message, verbosity_level=0):
        if verbosity_level <= self.max_verbosity:
            print(message)

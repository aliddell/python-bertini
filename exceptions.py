class BertiniError(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

class NoBertiniError(Exception):

    def __init__(self):
        self.message = "You don't seem to have Bertini installed anywhere " \
                       "I can find it."
    def __str__(self):
        return self.message
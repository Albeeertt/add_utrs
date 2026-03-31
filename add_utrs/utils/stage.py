

class Stage:

    def __init__(self, func, inputs, outputs):

        self.func = func
        self.inputs = inputs
        self.outputs = outputs

    def __call__(self, context):

        args = [context[key_out] for key_out in self.inputs]
        return_func = self.func(*args)

        if len(self.outputs) == 1:
            context[self.outputs] = return_func
        else:
            for key_out, only_return_func in zip(self.outputs, return_func):
                context[key_out] = only_return_func

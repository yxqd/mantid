

def _create_wrapper_function(name):
    """
    Create fake functions for the given name

    **It should not be called directly**

    :param name: name of fake function
    """
    # ------------------------------------------------------------------------------------------------
    def wrapper_function(*args, **kwargs):
        name_to_constructor = {
            'CompositeFunction': CompositeFunctionWrapper,
            'ProductFunction': ProductFunctionWrapper,
            'Convolution': ConvolutionWrapper,
            'MultiDomainFunction': MultiDomainFunctionWrapper,
            }
        # constructor is FunctionWrapper if the name is not in the registry.
        if name in name_to_constructor:
            return name_to_constructor[name](*args, **kwargs)
        return FunctionWrapper(name, **kwargs)

    # ------------------------------------------------------------------------------------------------
    wrapper_function.__name__ = name
    globals()[name] = wrapper_function

fnames = FunctionFactory.getFunctionNames()
for i, val in enumerate(fnames):
    _create_wrapper_function(val)

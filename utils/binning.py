def rebin(array, new_shape):

    '''
    
    '''

    # determine if 1D or 2D:
    arr_type = [1, 2]
    old_shape = len(array.shape)

    assert old_shape in arr_type

    if old_shape == 1: # 1D

        shape = array.shape[0] // new_shape[0]
        return array.reshape(-1, shape).mean(axis = 1)
   
    elif old_shape == 2: 

        shape = (new_shape[0], array.shape[0] // new_shape[0],
                new_shape[1], array.shape[1] // new_shape[1])
        
        return array.reshape(shape).mean(-1).mean(1)

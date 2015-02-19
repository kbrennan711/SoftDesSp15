"""Header Comment:
Kelly Brennan
Software Design
Professor Paul Ruvolo
February 18th, 2015

This code generates computational art by using recursion to build a random mathematical function. 
The combination of the evaluated functions for the red, green and blue channels produces the final image. 
A unique pictures is saved in the current directory each time the code is run. 
"""

import random
from math import sin, cos, pi, tan
from PIL import Image


#Random functions used


def build_random_function(min_depth, max_depth):
    """ Builds a random function of depth at least min_depth and depth
        at most max_depth (see assignment writeup for definition of depth
        in this context)

        min_depth: the minimum depth of the random function
        max_depth: the maximum depth of the random function
        returns: the randomly generated function represented as a nested list
                 (see assignment writeup for details on the representation of
                 these functions)
    """
    func_list = ['x', 'y', 'prod', 'avg', 'cos_pi', 'sin_pi', 'sqr', 'tan']
    if min_depth > 0:
        function = random.choice(func_list[2:])
    elif max_depth <= 0:
        function = random.choice(func_list[:2])
    else:
        function = random.choice(func_list)
    if function == 'x':
        f = lambda x,y: x
    elif function == 'y':
        f = lambda x,y: y
    elif function == 'prod': 
        nested_func = build_random_function(min_depth-1, max_depth-1)
        nested_func2 = build_random_function(min_depth-1, max_depth-1)
        f = lambda x,y: nested_func(x,y)*nested_func2(x,y)
    elif function == 'avg':
        nested_func = build_random_function(min_depth-1, max_depth-1)
        nested_func2 = build_random_function(min_depth-1, max_depth-1)
        f = lambda x,y: 0.5*(nested_func(x,y)*nested_func2(x,y))
    elif function == 'cos_pi':
        nested_func = build_random_function(min_depth-1, max_depth-1)
        f = lambda x,y: cos(pi*nested_func(x,y))
    elif function == 'sin_pi':
        nested_func = build_random_function(min_depth-1, max_depth-1)
        f = lambda x,y: sin(pi*nested_func(x,y))
    elif function == 'tan':
        nested_func = build_random_function(min_depth-1, max_depth-1)
        f = lambda x,y: tan(pi*nested_func(x,y))
    elif function == 'sqr':
        nested_func = build_random_function(min_depth-1, max_depth-1)
        f = lambda x,y: nested_func(x,y)**2
    return f


def remap_interval(val, input_interval_start, input_interval_end, output_interval_start, output_interval_end):
    """ Given an input value in the interval [input_interval_start,
        input_interval_end], return an output value scaled to fall within
        the output interval [output_interval_start, output_interval_end].

        val: the value to remap
        input_interval_start: the start of the interval that contains all
                              possible values for val
        input_interval_end: the end of the interval that contains all possible
                            values for val
        output_interval_start: the start of the interval that contains all
                               possible output values
        output_interval_end: the end of the interval that contains all possible
                            output values
        returns: the value remapped from the input to the output interval

        >>> remap_interval(0.5, 0, 1, 0, 10)
        5.0
        >>> remap_interval(5, 4, 6, 0, 2)
        1.0
        >>> remap_interval(5, 4, 6, 1, 2)
        1.5
    """
    a = float(input_interval_end - val)
    b = float(val - input_interval_start)
    c = b/(a+b)
    output = c*(output_interval_end - output_interval_start) + output_interval_start
    return output


def color_map(val):
    """ Maps input value between -1 and 1 to an integer 0-255, suitable for
        use as an RGB color code.

        val: value to remap, must be a float in the interval [-1, 1]
        returns: integer in the interval [0,255]

        >>> color_map(-1.0)
        0
        >>> color_map(1.0)
        255
        >>> color_map(0.0)
        127
        >>> color_map(0.5)
        191
    """
    # NOTE: This relies on remap_interval, which you must provide
    color_code = remap_interval(val, -1, 1, 0, 255)
    return int(color_code)


def generate_art(filename, x_size=350, y_size=350):
    """ Generate computational art and save as an image file.

        filename: string filename for image (should be .png)
        x_size, y_size: optional args to set image dimensions (default: 350)
    """
    # # Functions for red, green, and blue channels - where the magic happens!
    red_function = build_random_function(7,15)
    green_function = build_random_function(7,15)
    blue_function = build_random_function(7,15)

    # Create image and loop over all pixels
    im = Image.new("RGB", (x_size, y_size)) #Image size
    pixels = im.load()
    for i in range(x_size): #Nested loop where pixel value is set
        for j in range(y_size): #based on evaluating R, G, and B channel functions
            x = remap_interval(i, 0, x_size, -1, 1)
            y = remap_interval(j, 0, y_size, -1, 1)
            pixels[i, j] = (
                color_map(red_function(x,y)),
                color_map(green_function(x,y)),
                color_map(blue_function(x,y))
                )
    im.save(filename)

# # 
if __name__ == '__main__':
    import doctest
    # doctest.run_docstring_examples(remap_interval, globals())
    doctest.testmod()

    # Create some computational art!
    # TODO: Un-comment the generate_art function call after you
    #       implement remap_interval and evaluate_random_function
    #for i in range(15):
    generate_art("myart_part2_lambdaTest8.png")

    # Test that PIL is installed correctly
    # TODO: Comment or remove this function call after testing PIL install
    # test_image("noise.png")
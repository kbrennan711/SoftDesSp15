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
from math import sin, cos, pi
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
    #Scheme for representing function compositions ['elementary function name here', argument 1, argument 2]
    func_list = ['x', 'y', 'prod', 'avg', 'cos_pi', 'sin_pi']
    if min_depth > 0:
        f = random.choice(func_list[2:])
    elif max_depth <= 0:
        f = random.choice(func_list[:2])
        return [f] #Makes it a list instead of a string
    else:
        f = random.choice(func_list)
    if f == 'x':
        nested_func = lambda x,y: x(x,y)
    elif f == 'y':
        nested_func = lambda x,y: y(x,y)
    elif f == 'prod': 
        nested_func = (lambda x,y: x)(x,y)*(lambda x,y: y)(x,y)
    elif f == 'avg':
        nested_func = (lambda x,y: x)(x,y)+(lambda x,y: y)(x,y) #note have not found the mean... yet
    elif f == 'cos_pi':
        nested_func = lambda x,y: cos(pi*x)(x,y)
    elif f == 'sin_pi':
        nested_func = lambda x,y: sin(pi*x)(x,y)
    return [nested_func, build_random_function(min_depth-1, max_depth-1)]

# green = build_random_function(2,3)
# green_channel_pixel_for_x_y = green(1,1)
# print build_random_function(2,3)


def evaluate_random_function(f, x, y):
    """ Evaluate the random function f with inputs x,y
        Representation of the function f is defined in the assignment writeup

        f: the function to evaluate
        x: the value of x to be used to evaluate the function
        y: the value of y to be used to evaluate the function
        returns: the function value

        >>> evaluate_random_function(['x'],-0.5, 0.75)
        -0.5
        >>> evaluate_random_function(['y'],0.1,0.02)
        0.02
    """
    if f[0] == 'x': 
        return x
    elif f[0] == 'y':
        return y
    elif f[0] == 'prod':
        return evaluate_random_function(f[1], x, y)*evaluate_random_function(f[2], x, y)
    elif f[0] == 'avg':
        return  0.5*(evaluate_random_function(f[1], x, y)+evaluate_random_function(f[2], x, y))
    elif f[0] == 'cos_pi':
        return cos(pi*evaluate_random_function(f[1], x, y))
    elif f[0] == 'sin_pi':
        return sin(pi*evaluate_random_function(f[1], x, y))


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
    # red_function = build_random_function(7,15)
    # green_function = build_random_function(7,15)
    # blue_function = build_random_function(7,15)

    # print red_function
    # print green_function
    # print blue_function

    # Create image and loop over all pixels
    im = Image.new("RGB", (x_size, y_size)) #Image size
    pixels = im.load()
    for i in range(x_size): #Nested loop where pixel value is set
        for j in range(y_size): #based on evaluating R, G, and B channel functions
            x = remap_interval(i, 0, x_size, -1, 1)
            y = remap_interval(j, 0, y_size, -1, 1)
            pixels[i, j] = (
                    color_map(evaluate_random_function(red_function, x, y)),
                    color_map(evaluate_random_function(green_function, x, y)),
                    color_map(evaluate_random_function(blue_function, x, y))
                    ) #Obtain intensity for each color channel

    im.save(filename)


if __name__ == '__main__':
    import doctest
    # doctest.run_docstring_examples(remap_interval, globals())
    doctest.testmod()

    # Create some computational art!
    # TODO: Un-comment the generate_art function call after you
    #       implement remap_interval and evaluate_random_function
    #for i in range(15):
    # generate_art("myart_part2_lambdaTest.png")

    # Test that PIL is installed correctly
    # TODO: Comment or remove this function call after testing PIL install
    # test_image("noise.png")
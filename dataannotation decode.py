# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 17:23:12 2023

@author: mlaur
"""




def create_pyramid(nums):
    """
    Creates a pyramid structure based on a list of numbers.

    Args:
        nums (list): A list of numbers.

    Returns:
        list: A list of lists representing the pyramid structure.
              Each sublist contains numbers from the input list, forming a level of the pyramid.
              Returns False if unable to create a complete pyramid.
    """
    level_length = 0   # Initialize the level length variable
    pyramid = []     # Initialize the list to store the pyramid structure
    input_list = sorted(nums).copy()  # Create a sorted copy of the input list to avoid modifying the original

    while len(input_list) > 0:  # Continue until all numbers are used
        level_length += 1  # Increment the level length for each iteration
        if len(input_list) < level_length:
            return False  # Return False if there are not enough numbers for the current level 
        level = input_list[:level_length]  # Extract numbers for the current level 
        input_list = input_list[level_length:]  # Remove the used numbers from the input list

        pyramid.append(level)  # Append the current level to the pyramid structure

    return pyramid  # Return the completed pyramid structure


def decode(message_file):
    """
    Reads an encoded message from a .txt file and returns the decoded version as a string.

    Args:
        message_file (str): The path to the .txt file containing the encoded message.

    Returns:
        str: The decoded message as a string.
        False: If the pyramid creation is unsuccessful.
    """
    dictionary = {}  # Initialize a dictionary to store number-word associations

    # Read and parse the .txt file
    with open(message_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.split()
            dictionary[int(parts[0])] = parts[1]
    
    numbers = list(dictionary.keys())
    pyramid = create_pyramid(numbers)  # Create a pyramid based on the list of numbers
    
    # Return False if pyramid creation is unsuccessful
    if not pyramid:
        return False

    words = [dictionary[level[-1]] for level in pyramid]  # Extract words based on pyramid structure
    decoded_message = ' '.join(words)  # Join the words to form the decoded message
    
    return decoded_message  # Return the decoded message

# Example usage:
message_file = 'message_file.txt'  # Replace with the actual file path
result = decode(message_file)
print(result)

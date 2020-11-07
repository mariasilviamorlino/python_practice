"""
from https://www.w3resource.com/python-exercises/class-exercises/
'write a Python class to convert an integer to a roman numeral'

strategy (integer to roman): left-pad the given integer (<3999) so it has 4 digits
transform every digit into its roman number equivalent"""


def left_pad(num):
    num = str(num)
    pad = '0'*(4-len(num))
    num = pad+num
    return num


def digit_to_roman(num, magnitude):
    magnitude_mapping = [['I', 'V', 'X'],
                ['X', 'L', 'C'],
                ['C', 'D', 'M'],
                ['M', 'n', 'o']]   
    digit_mapping = ['', 'W', 'WW', 'WWW',
                    'WY', 'Y', 'YW',
                    'YWW', 'YWWW', 'WZ']
    
    digit = digit_mapping[num]
    magnitude_chars = magnitude_mapping[magnitude]
    digit = digit.replace('W', magnitude_chars[0])
    digit = digit.replace('Y', magnitude_chars[1])
    digit = digit.replace('Z', magnitude_chars[2])
    return digit


def int_to_roman(num):
    if num > 3999:
        print("sorry, can't transform this integer to roman!")
        return
    num = left_pad(num)
    roman = ''
    mag = 3
    for i in num:
        roman += digit_to_roman(int(i), mag)
        mag -= 1
    return roman


if __name__ == '__main__':
    print(int_to_roman(1))
    print(int_to_roman(3999))
    print(int_to_roman(27))


# notes (after peeking solutions)
# failed to think about modular arithmetics
# still think my implementation is not horrible

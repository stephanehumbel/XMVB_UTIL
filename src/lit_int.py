#!/usr/bin/env python3
import re
import routines
def detect_keyword(line, keyword):
    if keyword in line:
        return True
    return False

def read_integer_after_keyword(line, keyword):
    words = re.split(' |=|\n',line)
    print('1',words)
#    words = [word.split('=') for word in words]
#    words = [word.split('p') for word in words]
#    print('2',words)
    index = words.index(keyword)
    integer = int(words[index + 1])
    return integer

# Example usage
list= [1,2,3,4,6,7,9,8,10]
print(list)
list.sort()
print(list)
str1=routines.makeSTR(list)
print (str1)
quit()
filename=input("Enter the file name: ")
keyword = input("Enter a keyword: ")
with open(filename , "r") as file:
    for line in file:
        if detect_keyword(line, keyword):
            integer = read_integer_after_keyword(line, keyword)
            print("The integer after the keyword is:", integer)
            break
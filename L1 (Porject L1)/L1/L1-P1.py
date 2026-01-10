#   Laboratory 1 Problem 1
#   Make an application that is able to find the alphabet of a sequence of text. The sequence may be an ARN or ADN sequence or proton sequence.

def find_alphabet(text):
    letters_only=[]
    for i in text:
        if i.isalpha():
            letters_only+=[i.lower()]
    unique=[]
    for i in letters_only:
        if i not in unique:
            unique+=[i]
    unique.sort()
    a=''
    for i in unique:
        a+=i
    return a

def main():
    t=input("Enter a text sequence: ")
    print("The alphabet of the text is: "+find_alphabet(t))

main()
